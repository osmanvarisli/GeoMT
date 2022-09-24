"""
@author: Osman VARIŞLI
"""
#from bs4 import BeautifulSoup
import requests
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import io
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import Volcano as VLC # bizim class



from bioinfokit import  visuz

#from statsmodels.stats.proportion import proportions_ztest
from scipy.stats import ttest_ind
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

import mygene
#user_path="UsersFiles/osman/"

from flask import Flask,render_template,request,send_file,session
from logging import FileHandler,WARNING
app = Flask(__name__)

file_handler=FileHandler('errorlog.txt')
file_handler.setLevel(WARNING)
app.logger.addHandler(file_handler)

app.secret_key = "abc" 

def replaceNull(x,val='0'):
    if x is None: return val
    if x=='': return val
    else: return x

def dosya_bilgileri(acc_Code,Platform):
    """
    Platformalar olduğunda hata veriyor.
        https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1460/matrix/
        https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181735/matrix/

        yukarıda ikisini karşılaştır.
    """
    user_path=session['user_path']

    if Platform !='':
       Platform='-'+Platform 
    category=acc_Code[:3]
    if category=='GSE':#Series
        s_code=acc_Code[:len(acc_Code)-3]
        url="https://ftp.ncbi.nlm.nih.gov/geo/series/"+s_code+"nnn/"+acc_Code+"/matrix/"+acc_Code+Platform+"_series_matrix.txt.gz"
        file_name=user_path+acc_Code+Platform+"_series_matrix.txt"
    elif  category=='GDS':#DataSet
        s_code=acc_Code[:4]
        url="https://ftp.ncbi.nlm.nih.gov/geo/datasets/"+s_code+"nnn/"+acc_Code+"/soft/"+acc_Code+".soft.gz"
        file_name=user_path+acc_Code+".soft"
    else:
        return '','',''
    #print(url)
    return s_code,url,file_name
def csv_oku(dosya,delimiter):
    user_path=session['user_path']
    return pd.read_csv(user_path+dosya,delimiter=delimiter)

def csv_kaydet(df,dosya):
    user_path=session['user_path']
    return df.to_csv(user_path+dosya) 

def dosyaoku(dosya_yol):
    user_path=session['user_path']
    def dosya_ayikla(dosya_yol):
        d=dosya_yol.split('/')
        dosya=d[2]
        category=dosya[:3]
        if category=='GSE':#Series
            baslama_="series_matrix_table_begin"
            bitisi_="series_matrix_table_end"
        elif category=='GDS':#dataset
            baslama_="dataset_table_begin"
            bitisi_="dataset_table_end"

        yazmaya_basla=False
        with open(dosya_yol) as f:
            with open(user_path+"ham_", "w") as f1:
                for line in f:
                    if line.find(bitisi_)>0: yazmaya_basla=False
                    if yazmaya_basla: f1.write(line)
                    if line.find(baslama_)>0: yazmaya_basla=True

    
    dosya_ayikla(dosya_yol) 
    df=csv_oku("ham_",delimiter="\t")#dosyayı okumak için fonksiyorunu çağırmıyoruz.
    #df=csv_oku("ham_"+dosya,delimiter="\t")
    return df
    #return render_template('simple.html',  tables=[df.to_html(classes='data')], titles=df.columns.values)
def kayitlari_goster(grup,kac_satir,satir_no,query='',xsort='',id_ref='',xsort_asc_desc=''):
    #csv bul ve oku
    df=csv_oku(grup+".csv",delimiter=",")
    #sorgu varsa bul
    if query!='' :
        df=df.query(query) 
    #null değerleri 0 ile değiştir.
    df = df.fillna(0)

    #sıralama yapılmışsa kontrol et.
    if xsort is not None and xsort!='' and  xsort!='0':
        if xsort_asc_desc is not None and xsort_asc_desc!='' and  xsort_asc_desc!='0' and xsort_asc_desc=='ASC':
            df=df.sort_values(xsort, ascending=True)
            
        else:
            df=df.sort_values(xsort, ascending=False)
    print(xsort_asc_desc)
    if id_ref is not None and id_ref!='' and  id_ref!='0':
        #df=df.set_index('ID_REF').filter(regex=id_ref, axis=0)
        df=df[df['ID_REF'].str.contains(id_ref)]
    #ajax load dan gelmiş olabilir.
    df = df.iloc[satir_no:kac_satir+satir_no]
    return df

def gen_adlari(df):
    array_list=df['ID_REF'].values.tolist()
    mg = mygene.MyGeneInfo()
    df=mg.querymany(array_list, scopes='reporter',fields='query,ensembl.gene,symbol', as_dataframe=True)
    
    micro_RNA = ['SNO', 'ATP', 'MIR', 'LINC','PI','SI']
    x=df["symbol"].str.startswith(tuple(micro_RNA))
    df=df[x==False]


    if 'symbol' not in df.columns: df['symbol']='None'
    return df

def id_ref_gen_ad_birlestir(df):

    #df boş gelebilir kontrolu sağlanıyor.
    if len(df)<=0 :
        return df

    genler=gen_adlari(df)
    yeni_df=df.join(genler['symbol'], on='ID_REF')  
    
    yeni_df['symbol'] = yeni_df['symbol'].fillna('None')
    yeni_df['ID_REF']=yeni_df['ID_REF']+' ( '+ yeni_df['symbol']+' )'
 
    yeni_df=yeni_df.drop(['symbol'], axis=1)

    return yeni_df
def get_grp_json():
    user_path=session['user_path']
    f = open(user_path+"grp_json.txt", "r")
    grp_json=f.read() 
    return grp_json

def fold_change_query(column_names,down,up,max_p_value):
    #foldchange için sorgu oluştur
    qq="("
    for i in column_names:
        qq=qq+"(`"+i+"_P_value`<= "+max_p_value+") or "
    qq=qq[:-3] 
    qq=qq+") and ("
    for i in column_names:
        qq=qq+"(`"+i+"`<= "+down +" or `"+i+"` >= "+up+") or "

    qq=qq[:-3] 
    qq=qq+")"
    print(qq)   
    return qq
def permitasyon_cek(df,column_names):
    from itertools import permutations
    coll=column_names
    perm = permutations(coll[0:])

    perm_list=list(perm)
    q_list=[]

    for l in perm_list:
        q=''
        for j in l:
            q+='`'+j+'`>'
        q=q[0:len(q)-1]
        if not df.query(q).empty:
            q_list.append(q)
    return q_list   

def normalizasyon_veri(acc_Code,Platform,grp_json):
    try: 
        s_code,url,file_name=dosya_bilgileri(acc_Code,Platform)
    except:
        #hata oldu
        pass
    #Her seferinde grp_json post etmeme gerek kalmayacak.Dosyadan çekecem
    user_path=session['user_path']
    f = open(user_path+"grp_json.txt", "w")
    f.write(grp_json)
    f.close()

    datastore = json.loads(grp_json)
    i=0
    grupAd_list=[]
    d_df=[]
    d_df_sutun_toplam=[]
    d_df_genel_toplam=[]
    d_df_column_names=[]
    d_df_row_data=[]
    df_ortalamalar = pd.DataFrame()  
    
    butun_stutun_toplam=0
    sutunlar_say=0
    for i in range(len(datastore["pj"])):
        sutunlar=datastore["pj"][i]["grup"].split(' ,')
        sutunlar_say=sutunlar_say+len(sutunlar)
        aha=dosyaoku(file_name)[sutunlar].round(5)

        d_df_sutun_toplam.append(aha.sum())
        butun_stutun_toplam=butun_stutun_toplam+aha.sum().sum()
    print(butun_stutun_toplam)
    ortalama=butun_stutun_toplam/sutunlar_say
    
    renklencek_alanlar=[]
    for i in range(len(datastore["pj"])):
        sutunlar=datastore["pj"][i]["grup"].split(' ,')
        sutunlar.insert(0,"ID_REF")

        grupAd=datastore["pj"][i]["grupAd"]
        grupAd_list.append(grupAd)
        pandas_data=dosyaoku(file_name)[sutunlar].round(5)
        id_refler=pandas_data['ID_REF']
        sut_say=len(sutunlar)
        say=0

        yeni_alanlar=[]
        for j in reversed(sutunlar):
            say+=1
            if j!="ID_REF":
                #pandas_data[j+'__'] = pandas_data[j]
                pandas_data.insert(sut_say-say+1, 'Nom_'+str(j), (pandas_data[j]*(ortalama/d_df_sutun_toplam[i][j])).round(5))
                yeni_alanlar.append('Nom_'+str(j))
                renklencek_alanlar.append('Nom_'+str(j))
        
        pandas_data['C']=pandas_data[yeni_alanlar].mean(axis=1).round(5)
        df_ortalamalar[grupAd]=pandas_data['C']
        csv_kaydet(pandas_data,grupAd+'.csv')


        d_df.append(pandas_data)
        d_df_column_names.append(d_df[i].columns.values)
        d_df_row_data.append(list(d_df[i].values.tolist()))

    # ortalama sekmesi için verileri ekliyorum
    grupAd_list.append('Averages')
    df_ortalamalar.insert(0,'ID_REF',id_refler)
    d_df_column_names.append(df_ortalamalar.columns.values)
    d_df_row_data.append(list(df_ortalamalar.values.tolist()))
    csv_kaydet(df_ortalamalar,'normalize_ortalamalar.csv')
    
    return  grupAd_list,d_df_column_names,d_df_row_data,d_df_sutun_toplam,butun_stutun_toplam,ortalama,renklencek_alanlar

def basit_normalizasyon(df):

    def sutun_bul(df):
        sutunlar=[]    
        for i in df.columns.values:
            olmamali=['Unnamed: 0','ID_REF','IDENTIFIER']
            if i not in olmamali:
                sutunlar.append(i)
        return sutunlar
    s=sutun_bul(df)
    toplamlar=df[s].sum()
    print(toplamlar.sum())
    print('osfadsfadsf')
    ortalama=toplamlar.sum()/len(df[s].columns)
    for i in s:
        normalizasyon_katsayisi=ortalama/df[i].sum()
        df[i]=df[i]*normalizasyon_katsayisi
        df = df.rename({i: 'Norm_'+i}, axis=1)
    return df

@app.route('/', methods=['POST','GET'])
def ana_sayfa():
    def temp_sil():
        import datetime

        path = r'C:\Users\PC\Desktop\Flask\UsersFiles'

        today = datetime.datetime.today()
        os.chdir(path) 
        print(os.getcwd())
        for root,directories,files in os.walk(path,topdown=False): 
            print(files)
            for name in files:
                print(name)
                t = os.stat(os.path.join(root, name))[8] 
                filetime = datetime.datetime.fromtimestamp(t) - today
                print(os.path.join(root, name), filetime.days)    
                if filetime.days <= -1:
                    print(os.path.join(root, name), filetime.days)
                    os.remove(os.path.join(root, name))

    if session.get('user_path') is None:
        session['user_path']= "UsersFiles/"+str(np.random.randint(9999999))+'/'

    user_path=session['user_path']
    if not os.path.exists(user_path):
        os.mkdir(user_path)
    
    #temp_sil()# yedi günden fazla olan usersfileleri sil
        
    return render_template('index.html')

@app.route('/data_list/', methods=['POST','GET'])
def ham_data_list():
    def dosyaindir(url,path):
        filename = path+url.split("/")[-1]
        with open(filename , 'wb') as f:
            r = requests.get(url)
            f.write(r.content)
    def unzip(dosya):
        import gzip
        import shutil
        with gzip.open(dosya, 'rb') as f_in:
            with open(dosya[:len(dosya)-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    user_path=session['user_path']
    if request.method=='POST':
        acc_Code=request.form['acc_Code'].strip()
        Platform=request.form['Platform'].strip()
        s_code,url,file_name=dosya_bilgileri(acc_Code,Platform)
        try:
            err='File Not Downloaded.'
            dosyaindir(url,user_path)
            err='Zipped file could not be opened'
            unzip(file_name+".gz")
            #zipten çıktıktan sonra dosyayı silebiliriz.

            err='The file could not be read, there may be no data in it'
            df=dosyaoku(file_name).head(100)
            
            #bütün listeyi eklememiz için baştan okuyoruz.

            err=''
        except:
            #raise
            return render_template('error.html',err=err)
        else:  
            parametreler={
                'acc_Code':acc_Code,
                'Platform':Platform,
                'column_names':df.columns.values,
                'row_data':list(df.round(5).values.tolist()),
            }   
            return  render_template('ham_data_list.html', parametreler=parametreler,
                           zip=zip)
    return render_template('error.html',err='Wrong page')
@app.route('/normalization/', methods=['POST','GET'])
def normalization():
    user_path=session['user_path']
    acc_Code=request.form['acc_Code'].strip()
    Platform=request.form['Platform'].strip()
    s_code,url,file_name=dosya_bilgileri(acc_Code,Platform)
    df=dosyaoku(file_name)
    df_norm=basit_normalizasyon(df)
    parametreler={
        'acc_Code':acc_Code,
        'Platform':Platform,
        'column_names':df.columns.values,
        'norm_column_names': df_norm.columns.values ,
        'norm_row_data':list(df_norm.head(100).round(5).values.tolist()),
    }   
    return  render_template('simple_normalization.html', parametreler=parametreler,
                    zip=zip)   

@app.route('/group/', methods=['POST','GET'])
def group():
        
        acc_Code=request.form['acc_Code']
        grp_json=request.form['grp_json']
        Platform=request.form['Platform']
        grupAd_list,d_df_column_names,d_df_row_data,d_df_sutun_toplam,butun_stutun_toplam,ortalama,renklencek_alanlar=normalizasyon_veri(acc_Code,Platform,grp_json)

        return  render_template('group.html', grup_list=grupAd_list,column_names=d_df_column_names, row_data=d_df_row_data,acc_Code=acc_Code,grp_json=grp_json,
                            renklencek_alanlar=renklencek_alanlar,zip=zip)

@app.route('/ajax_load/', methods=['POST','GET'])
def ajax_load():
    #acc_Code=request.form['acc_Code']
    #grp_json=request.form['grp_json']
    grp_json=get_grp_json()
    satir_no=request.form['no']
    
    grup_sec=int(request.form['grup_sec'])
    
    grp_json = json.loads(grp_json)
    #datastore = json.loads(grp_json)

    if grup_sec>(len(grp_json["pj"])-1):#gruplandırılmış normalizasyondan gelip gelmediğini kontrol edecez.  
        grup_ad='normalize_ortalamalar'
    else: 
        grup_ad=grp_json["pj"][grup_sec]["grupAd"]
    #print("************************************")
    #print (str(grup_ad))
    #print("************************************")
    df=kayitlari_goster(grup_ad,500,int(satir_no)).round(5)
    
    row_data=list(df.values.tolist())
    #print(df)

    return  render_template('table_ajax_load.html',row_data=row_data)

@app.route('/fold_change_ajax_load/', methods=['POST','GET'])
def fold_change_ajax_load():
    #fold change için scroll en sonra geldiğinde ve ilk yükleme için kullanıldı.
    print(request.form)
    grup_ad=request.form['grup_ad']
    satir_no=request.form['no']
    down=request.form['down']
    up=request.form['up']
    max_p_value=request.form['max_p_value']

    if request.form['xsort'] is None: xsort=''
    else: xsort=request.form['xsort']

    if request.form['xsort_asc_desc'] is None: xsort_asc_desc=''
    else: xsort_asc_desc=request.form['xsort_asc_desc']

    if request.form['id_ref'] is None:
        id_ref=''
    else:
        id_ref=request.form['id_ref']

    gelen_yer=replaceNull(request.form.get('gelen_yer'))#popup penceresinden gelmiş olabilir

    column_names=request.form['column_names']
    column_names=column_names.replace(" ",",")
    column_names=column_names.replace(",,",",")
    column_names= eval(column_names.split(' ')[0])

    qq=fold_change_query(column_names,down,up,max_p_value)
    if gelen_yer=='alt_sorgu':
        sorgu=replaceNull(request.form.get('sorgu'))
        qq='( '+qq+' ) and '+sorgu
    else:
        sorgu=''
    #print(qq)   
    df=kayitlari_goster(grup_ad,500,int(satir_no),qq,xsort,id_ref,xsort_asc_desc).round(5)
    
    sutunlar=getir_fold_change_sutunlar(df)#sadece foldchange sutunları
    say=0
    sutun_index=[]
    for val in df.columns.values:
        say=say+1
        if val not in sutunlar:
            sutun_index.append(say) #foldchange dışında ki sutunları bul

    df=id_ref_gen_ad_birlestir(df)
    #print(request.form['grp_json'])
    row_data=list(df.values.tolist())

    return render_template('fold_change_table_ajax_load.html',sutun_index=sutun_index,row_data=row_data,sira=4,down=down,up=up,sorgu=sorgu)

def getir_fold_change_sutunlar(df):
    #
    #Fold_change tablosu içinde sadece ana kolonları döndürür. P_value,ID_REF gibi sutnları döndürmez.
    #
    sutunlar=[]    
    for i in df.columns.values:
        if i[len(i)-7:]!='P_value':
            olmamali=['Unnamed: 0','ID_REF']
            if i not in olmamali:
                sutunlar.append(i)
    return sutunlar

@app.route('/fold_change/', methods=['POST','GET'])
def fold_change():
    down=str(float(replaceNull(request.form.get('down'),'-2')))
    up=str(float(replaceNull(request.form.get('up'),'2')))
    max_p_value=str(float(replaceNull(request.form.get('max_p_value'),'0.00005')))
    column_names=replaceNull(request.form.get('column_names'))
    xsort=replaceNull(request.form.get('xsort'))
    xsort_asc_desc=replaceNull(request.form.get('xsort_asc_desc'))

    if request.form.get('id_ref') is None:
        id_ref=''
    else:
        id_ref=request.form.get('id_ref')

    control=request.form['control']
    #grp_json=request.form['grp_json']
    grp_json=get_grp_json()
    gelen_yer=replaceNull(request.form.get('gelen_yer'))#popup penceresinden gelmiş olabilir

    datastore = json.loads(grp_json)
    df_control=csv_oku(control+".csv",delimiter=",")
    
    #control alanın sutunların adlarını çekiyoruz. 
    for i in range(len(datastore["pj"])):
        grup=datastore["pj"][i]["grupAd"]
        if (grup==control):
            control_sutunlar=datastore["pj"][i]["grup"].split(' ,')
    control_sutunlar = ['Nom_'+x for x in control_sutunlar]
    
    if gelen_yer=='alt_sorgu' or gelen_yer=='filter':
        df=csv_oku("foldchange.csv",delimiter=",")
        column_names=getir_fold_change_sutunlar(df)
        column_names_title=df.columns.values
        column_names_title = np.delete(column_names_title, 0) 
    else:
        df = pd.DataFrame()
        df=df.round(5)#virgülden sonra 5 değer
        df['ID_REF']=df_control['ID_REF']
        fold_change_sutunlar=[]
        for i in range(len(datastore["pj"])):
            grup=datastore["pj"][i]["grupAd"]
            csv_sutunlar=datastore["pj"][i]["grup"].split(' ,')
            csv_sutunlar = ['Nom_'+x for x in csv_sutunlar]
            #csv_sutunlar = [x for x in csv_sutunlar]

            df_csv=csv_oku(grup+".csv",delimiter=",") 

            if (grup!=control):
                fold_change_sutun=grup+'__'+control
                df[fold_change_sutun]=df_csv['C']/df_control['C']
                fold_change_sutunlar.append(fold_change_sutun)
                
                df = df.fillna(0)
                
                df_csv = df_csv.fillna(0)
                df_control = df_control.fillna(0)

                df_control['birlesik']=df_control[control_sutunlar].apply(lambda row: ';'.join(row.values.astype(str)), axis=1) 
                df_csv['birlesik']=df_csv[csv_sutunlar].apply(lambda row2: ';'.join(row2.values.astype(str)), axis=1)


                df['t']=df_csv['birlesik']+'_'+df_control['birlesik']

                f = lambda x: ttest_ind(np.asarray(x.split('_')[0].split(';'), dtype=np.float32),np.asarray(x.split('_')[1].split(';'), dtype=np.float32)).pvalue 
                df[grup+'__'+control+'_P_value'] = df['t'].apply(f)
                
                df = df.drop('t', 1)
                
        csv_kaydet(df,'foldchange.csv')
        column_names=fold_change_sutunlar
        column_names_title=df.columns.values
    
    if len(column_names)<2:#ID_REF de var 
        q_list=''
    else:
        qq=fold_change_query(column_names,down,up,max_p_value)
        q_list=permitasyon_cek(df.query(qq),column_names)
    


    if gelen_yer=='alt_sorgu':
        sorgu=replaceNull(request.form.get('sorgu'))
        return  render_template('fold_change_alt_sorgu.html',df=df,column_names=column_names,column_names_title=column_names_title,
                            down=down,up=up,control=control,xsort=xsort,xsort_asc_desc=xsort_asc_desc,id_ref=id_ref,max_p_value=max_p_value,
                            grp_json=grp_json,q_list=q_list,sorgu=sorgu,
                            zip=zip)
    else:
        return  render_template('fold_change.html',df=df,column_names=column_names,column_names_title=column_names_title,
                            down=down,up=up,control=control,xsort=xsort,xsort_asc_desc=xsort_asc_desc,id_ref=id_ref,max_p_value=max_p_value,
                            grp_json=grp_json,q_list=q_list,
                            zip=zip)

@app.route('/k_means_ajax_load/', methods=['POST','GET'])
def k_means_ajax_load():
    satir_no=int(request.form['no'])
    sorgu=request.form['sorgu']
    #down=request.form['down']
    #up=request.form['up']
    kac_satir=500
    df=csv_oku("K-means.csv",delimiter=",").round(5)
    df=df[(df.cluster == int(sorgu))]

    df = df.fillna(0)
    df = df.iloc[satir_no:kac_satir+satir_no]
    df=id_ref_gen_ad_birlestir(df)
    row_data=list(df.values.tolist())
    return  render_template('table_ajax_load.html',row_data=row_data)

@app.route('/k_means_group_table/', methods=['POST','GET'])
def k_means_group_table():
    
    sorgu=request.form['sorgu']
    down=str(float(replaceNull(request.form.get('down'))))
    up=str(float(replaceNull(request.form.get('up'))))    

    df=csv_oku("K-means.csv",delimiter=",")
    column_names=getir_fold_change_sutunlar(df)
    column_names_title=df.columns.values
    column_names_title = np.delete(column_names_title, 0) 

    parametreler={
        'column_names_title':column_names_title,
        'sorgu':sorgu,
        'down':down,
        'up':up,
    }   
    return  render_template('k_means_alt_sorgu.html', parametreler=parametreler,
                    zip=zip)

@app.route('/kmeans/', methods=['POST','GET'])
def kmeans():
    k_means_group_count=request.form['k_means_group_count']
    try:
        if int(k_means_group_count)<2 : k_means_group_count=2
    except:
        k_means_group_count=2

    df=csv_oku("foldchange.csv",delimiter=",")

    #fold change filtrele
    down=replaceNull(request.form.get('down'))
    up=replaceNull(request.form.get('up'))
    xsort=replaceNull(request.form.get('xsort'))
    max_p_value=replaceNull(request.form.get('max_p_value'))

    sutunlar=getir_fold_change_sutunlar(df)
    q=fold_change_query(sutunlar,down,up,max_p_value)

    df = df.drop('Unnamed: 0', 1)
    df=df.query(q) 
    # filtreleme bitti 


    df=df[sutunlar]
    df['Control'] = 1
    sutunlar.append('Control')
    #print(df)
    df = df.astype(np.float16)
    x=df.values.tolist()

       
    kmeans = KMeans(n_clusters=int(k_means_group_count))
    kmeans.fit(x)
    identified_clusters = kmeans.fit_predict(x)

    #print(identified_clusters)
    df['cluster'] = identified_clusters

    user_path=session['user_path']
    
    print(df)
    print('osman')
    plt.figure()
    pd.plotting.parallel_coordinates(df, 'cluster')
    plt.savefig(user_path+'/k-means.png')
    print('osman..')
    """     
    !!!!!!------Bunu Silme -----!!!
    Ayrıştırılmış K-means graafiği
    

    plt.figure()
    for i in range(int(k_means_group_count)):
        plt.subplot(2, 2,i+1)
        pd.plotting.parallel_coordinates(df[(df.cluster == i) ], 'cluster')

    plt.savefig(user_path+'/k-means2.png')
    """
    df_fold_change=csv_oku("foldchange.csv",delimiter=",")
    df_fold_change=df_fold_change.rename(columns={"Unnamed: 0": "No"})
    

    df=df.reset_index()
    df=df.rename(columns={"index": "No"})
    
    df=df[['No','cluster','Control']]

    column_names_title=df.columns.values

    df=df.set_index('No')
    df=df.join(df_fold_change, on='No')  
    
    df=df.drop(['No', 'Control'], axis=1)

    csv_kaydet(df,'K-means.csv')

    import random 
    parametreler={
        'k_means_group_count':int(k_means_group_count),
        'rnd':random.randint(0, 9999999999),
        'user_path':user_path,
    }   
    return  render_template('k-means.html', parametreler=parametreler,
                    zip=zip)
    
    #return str(df)
@app.route('/k-means.png', methods=['POST','GET'])
def k_means_plot():
    user_path=session['user_path']
    return send_file(user_path+'k-means.png', mimetype='image/png')

@app.route('/plot.png', methods=['POST','GET'])
def plot():
    
    if str(request.args.get('t'))=='satir':
        fig = satir_bazli_grafik(str(request.args.get('p'))) 
    elif str(request.args.get('t'))=='volcano':
        down=replaceNull(request.args.get('down'))
        up=replaceNull(request.args.get('up'))
        p=str(request.args.get('p'))
        fig = volcano_grafik(p,down,up)
    elif str(request.args.get('t'))=='normalization-kmeans':
        k_means_group_count=replaceNull(request.args.get('p'))
        """
        down=replaceNull(request.args.get('down'))
        up=replaceNull(request.args.get('up'))
        p=str(request.args.get('p'))
        """


        fig = normalization_kmeans_grafik(k_means_group_count)       
    elif str(request.args.get('t'))=='boxplot_foldchange':
        fig = boxplot_foldchange() 
    else:
        graph=str(request.args.get('g'))
        colons=str(request.args.get('p'))
        fig = create_figure(graph,colons)

    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')
def boxplot_foldchange():
    df=csv_oku("foldchange.csv",delimiter=",").columns[2:]

    fig, ax = plt.subplots()
    ax.set_title('Hide Outlier Points')
    ax.boxplot(df, showfliers=False)

    return fig

def normalization_kmeans_grafik(k_means_group_count):

    grp_json=get_grp_json()
    datastore = json.loads(grp_json)
    
    fig = Figure()
    
    #sutunların adlarını çekiyoruz. 
    say=0
    for i in range(len(datastore["pj"])):
        say=say+1
        grup=datastore["pj"][i]["grupAd"]
        sutunlar=datastore["pj"][i]["grup"].split(' ,')
        sutunlar = ['Nom_'+x for x in sutunlar]

        df=csv_oku(grup+".csv",delimiter=",")[sutunlar]
        df = df.fillna(0)

        pca = PCA(2)
        df = pca.fit_transform(df)
        kmeans = KMeans(n_clusters= int(k_means_group_count)) 
        label = kmeans.fit_predict(df)
        u_labels = np.unique(label)

        axis = fig.add_subplot(2, 2, say)
        axis.title.set_text(grup)
        for i in u_labels:
            axis.scatter(df[label == i , 0] , df[label == i , 1] , label = i)
        axis.legend()
        #axis.show()       

    return fig




def volcano_grafik(grup,down,up,max_p_value):  
    df=csv_oku("foldchange.csv",delimiter=",")

    qq=fold_change_query([grup],down,up,max_p_value)
    df=df.query(qq)
    
    df['log2']=np.log2(df[grup])
    df=df.dropna(subset=['log2'])

    volc = VLC.Volcano(df['log2'], df[grup+"_P_value"], df["ID_REF"], s_curve_x_axis_overplot=0.05, s_curve_y_axis_overplot=0.05)
    fig = volc.get_fig()
    return fig

def satir_bazli_grafik(satir):
    """
    user_path=session['user_path']
    f = open(user_path+"grp_json.txt", "r")
    grp_json=f.read() 

    datastore = json.loads(grp_json)
    
    fig = Figure()
    axis = fig.add_subplot(2, 1, 1)

    renk=[]
    y=[]
    x=[]
    r=['blue','green','red','orange','yellow','purple','pink']
    say=0
    for i in range(len(datastore["pj"])):
        say=say+1
        grup=datastore["pj"][i]["grupAd"]
        sutunlar=datastore["pj"][i]["grup"].split(' ,')
        sutunlar = ['Nom_'+x for x in sutunlar]
        sutunlar.append('ID_REF')
        for j in sutunlar : renk.append(r[say-1])

        df=csv_oku(grup+".csv",delimiter=",")[sutunlar]
        df = df.fillna(0)
        df=df.loc[df['ID_REF'] == satir]
        df=df.drop(columns=['ID_REF'])
        #print(df)
        for col in df.columns:
            #print('osman')
            #print(df.iloc[0][col])
            x.append(col)
            y.append(df.iloc[0][col]) 

    axis.tick_params(labelrotation=90)
    axis.bar(x, y, color = renk)
    return fig
    """
    df=csv_oku("foldchange.csv",delimiter=",")
    #df=df.loc[:,df.columns.str.endswith("P_value")]
    df=df.loc[df['ID_REF'] == satir]
    sutunlar=getir_fold_change_sutunlar(df)

    fig = Figure()
    axis = fig.add_subplot(1, 1, 1)
    xs = sutunlar#df.columns.values[2:]  

    ys = df[sutunlar].values[0]

    xs=list(xs)
    xs.append('Control')
    ys=list(ys)
    ys.append(1)
    axis.plot(xs, ys)
    return fig
   
def create_figure(graph,colons):
    #df=pd.read_csv("foldchange.csv",delimiter=",")
    df=csv_oku("foldchange.csv",delimiter=",")

    fig = Figure()
    axis = fig.add_subplot(1, 1, 1)
    """
    if graph=="grafik":
        #df_g=df.query(colons)
        xs = df[colon_arr[0]]
        ys = df[colon_arr[1]]       
        axis.plot(xs, ys)
    """

    return fig

def volcano_bioinfokit(grup,down,up):  
    df=csv_oku("foldchange.csv",delimiter=",")
    """
    qq=fold_change_query([grup],down,up)
    df=df.query(qq)
    print(qq)
    """
    df['log2']=np.log2(df[grup])
    df=df.dropna(subset=['log2'])
    user_path=session['user_path']
    visuz.gene_exp.volcano(df=df, lfc='log2', pv=grup+"_P_value",figname=user_path+'/Volcano_'+grup,figtype='png',sign_line=True)

@app.route('/volcano/<grup>/<rnd>/v.png', methods=['POST','GET'])
def vv(grup,rnd):
    user_path=session['user_path']
    return send_file(user_path+'Volcano_'+grup+'.png', mimetype='image/png')

@app.route('/grafik/volcano/<grup>', methods=['POST','GET'])
def volcano(grup):  
    down=replaceNull(request.form.get('down'))
    up=replaceNull(request.form.get('up'))

    volcano_bioinfokit(grup,down,up)
    return '<img src="../volcano/'+grup+'/'+str(np.random.randint(9999999))+'/v.png" width="100%">'
    #return '<img src="../plot.png?t=volcano&p='+grup+'&down='+down+'&up='+up+'" width="100%">'

@app.route('/grafik/normalization-kmeans/<grup>', methods=['POST','GET'])
def kmeans_g(grup):
    down=replaceNull(request.form.get('down'))
    up=replaceNull(request.form.get('up'))

    return '<img src="../plot.png?t=normalization-kmeans&p='+grup+'&down='+down+'&up='+up+'" width="100%">'

@app.route('/grafik/satir/<sorgu>', methods=['POST','GET'])
def grafik(sorgu):
    return '<img src="../plot.png?g=grafik&t=satir&p='+sorgu+'" width="100%">'

@app.route('/grafik/boxplot/', methods=['POST','GET'])
def boxplot():
    return '<img src="../plot.png?g=boxplot&t=boxplot_foldchange" width="100%">'    

@app.route('/export_csv/', methods=['POST','GET'])
def export_csv():
    
    dosya_ad=request.form['dosya_ad']
    df=csv_oku(dosya_ad+".csv",delimiter=",")

    if dosya_ad=='foldchange': 
        down=replaceNull(request.form.get('down'))
        up=replaceNull(request.form.get('up'))
        xsort=replaceNull(request.form.get('xsort'))
        
        max_p_value=replaceNull(request.form.get('max_p_value'))

        sutunlar=getir_fold_change_sutunlar(df)
        q=fold_change_query(sutunlar,down,up,max_p_value)

        df = df.drop('Unnamed: 0', 1)
        df=df.query(q) 
        if len(df)>0 :

            genler=gen_adlari(df)
            
            yeni_df=df.join(genler['symbol'], on='ID_REF') 
            df=yeni_df
            """
            yeni_df['symbol'] = yeni_df['symbol'].fillna('None')
            """
        #df=df.sort_values(xsort)
    elif dosya_ad=='k-means':
        sorgu=replaceNull(request.form.get('sorgu'))
        df=df[(df.cluster == int(sorgu))]

    csv_kaydet(df,'Download.csv')
    
    user_path=session['user_path']
    return send_file(user_path+'Download.csv', as_attachment=True)

@app.route('/list/', methods=['POST','GET'])
def ff():
    return 'loooo'

if __name__ == '__main__':
    app.debug = True
    app.run()
