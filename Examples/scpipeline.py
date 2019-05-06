import os, sys, csv,json,datetime,time,math,scipy.stats,collections,re;
from sklearn import preprocessing;
import numpy as np;
import pandas as pd;
import subprocess, gzip;
from subprocess import Popen, PIPE,STDOUT
import os.path;
import pymongo;
from bson.objectid import ObjectId;
from pymongo import MongoClient;
import scanpy;
import scanpy.api as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80)
from pymongo import IndexModel, ASCENDING, DESCENDING
import subprocess, gzip;
from subprocess import Popen, PIPE,STDOUT
import os.path;
import urllib.request

class ProcessPipline:
    def __init__(self):
        self.data="";
        self.CountsFile="";
        self.samples="";
        
    def runCellRange(self,workspace,fastq_path,samples,expect_cells,transcriptpath,run_name):
        run_name=run_name.replace(" ","_");
        resultfile ="";
        pcheck = subprocess.Popen("cellranger", shell=True, stdin=PIPE, stdout=PIPE,stderr=STDOUT)
        output =pcheck.stdout.read();
        output=str(output)
        if "command not found" in output:
            print("Please install cellranger 3.0");
            print("Tutorial: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation")
            print("If you have counts data. please skip this step.")
            return ""
        if not os.path.isdir(workspace):
            print("workspace is not a dir");
            return;

        if len(samples) ==0 :
            print("Please input samples");
            return;
        if not os.path.isdir(fastq_path):
            print("fastq path is not a dir");
            return;
        ResPath = workspace+"/"+run_name;
        if os.path.isdir(ResPath):
            print("run name is already in workspace");
            return;
        else:
            os.mkdir(ResPath);
        
        self.samples=samples;
        shstr = "";
        shstr +="wd=\"" +ResPath+  "\"\n"
        shstr +="cd ${wd}\n";
        
        for i in samples:
            jobstr="cellranger count ";
            jobstr+="--id "+i+" ";
            jobstr+="--fastqs "+fastq_path+" "
            jobstr+="--transcriptome "+transcriptpath+" ";
            jobstr+="--expect-cells "+ str(expect_cells)+" "
            jobstr+="--sample=\""+i+"\"";

            shstr+=jobstr+"\n";


        if len(samples) >1:
            csvstr=["library_id,molecule_h5"]
            for i in samples:
                temp=i+","+ResPath+"/"+i+"/outs/molecule_info.h5"
                csvstr.append(temp)

            csvstr="\n".join(csvstr);
            
            csvf=ResPath+"/"+run_name+".csv"
            with open(csvf,"w") as f:
                f.write(csvstr)
                
            shstr+="cellranger aggr "
            shstr+="--id aggr " 
            #--csv=test.csv --normalize=mapped
            shstr+= "--csv="+csvf+" --normalize=mapped"

            resultfile = "aggr";

        else:
            resultfile =samples[0]
            
        shfile=ResPath+"/"+run_name+".sh"
        with open(shfile,"w") as f:
            f.write(shstr)
            
        command = "bash "+shfile;
        
        print("it will take a few hours . please wait.....")
        prun = subprocess.Popen(command, shell=True, stdin=PIPE, stdout=PIPE,stderr=STDOUT)
        output =prun.stdout.read();
        print(output)
        print("---------------------------------------------------------------------------------");
        print("finish");
        resultpath = ResPath+"/"+resultfile;
        print("results path: "+resultpath);
        self.CountsFile=resultpath+"/outs/filtered_feature_bc_matrix/";
        
        print("counts file path: "+self.CountsFile);
        
        return self.CountsFile;
    
    def downloadTestData(self):
    
        datafolder="pbmc3k";
        downloadFile="bmc3k.tar.gz";
        dataname=datafolder+"/"+downloadFile;
        if os.path.exists(datafolder):
            pass;
        else:
            os.mkdir(datafolder);

        if os.path.exists(dataname):
            pass;
        else:
            urllib.request.urlretrieve("http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",dataname);
            print("download success");

        command = subprocess.Popen("cd "+datafolder+"\ntar -zxvf "+downloadFile, shell=True, stdin=PIPE, stdout=PIPE,stderr=STDOUT)
        output =command.stdout.read();
        countsFile=datafolder+"/filtered_gene_bc_matrices/hg19"
        self.CountsFile = countsFile;
        print(countsFile)
        return countsFile;


    def readData(self,countsFile=""):
        if countsFile=="":
            countsFile = self.CountsFile;
            
        if countsFile=="":
            print("please input counts file path");
            return ""
        
        self.CountsFile=countsFile;
        
        datapath = self.CountsFile;
        if os.path.isdir(datapath):
            files = os.listdir(datapath)
            for i in files:
                if i.endswith(".gz"):
                    print(i)
                    target = datapath+"/*.gz";
                    print(target)
                    command = subprocess.Popen("gunzip "+target, shell=True, stdin=PIPE, stdout=PIPE,stderr=STDOUT)
                    output =command.stdout.read();
                    break;
                    
            files=os.listdir(datapath);
            for i in files:
                if i =="features.tsv":
                    os.rename(datapath+"/features.tsv",datapath+"/genes.tsv");
                    break;
            files = list(os.listdir(datapath));
            if ('barcodes.tsv' in files) and ('barcodes.tsv' in files) and ("genes.tsv" in files):
                adata = sc.read_10x_mtx(datapath, var_names='gene_symbols');
                self.data=adata;
                self.preprocess();
            else:
                print("input data is not correct")
                return ""
            
        elif os.path.isfile(datapath):
            if datapath.endswith(".h5ad"):
                adata=sc.read(datapath);
            else:
                adata = sc.read_csv(datapath)
                adata = adata.T;
            self.data=adata;
            self.preprocess();
        else:
            print("file or dir not exists")
            return ""
        
    
    def preprocess(self):
        self.data.var_names_make_unique()
        sc.pp.filter_cells(self.data, min_genes=0);
        sc.pp.filter_genes(self.data, min_cells=1);
        mito_genes = [name for name in self.data.var_names if name.startswith('MT-')]
        self.data.obs['percent_mito'] = np.sum(self.data[:, mito_genes].X, axis=1)  / np.sum(self.data.X, axis=1)
        self.data.obs['n_counts'] = self.data.X.sum(axis=1);
        
        
    def QC(self,max_n_genes="" ,min_n_genes="",min_n_cells="",max_percent_mito=""):
        if min_n_genes!="":
            print("filter cells");
            sc.pp.filter_cells(self.data, min_genes=min_n_genes);
        if min_n_cells != "":
            print("filter genes");
            sc.pp.filter_genes(self.data, min_n_cells)
        if max_n_genes !="":
            print("filter n_genes < "+str(max_n_genes))
            self.data = self.data[self.data.obs['n_genes'] < max_n_genes , :]
        if max_percent_mito != "":
            print("filter percent_mito < "+str(max_percent_mito))
            self.data = self.data[self.data.obs['percent_mito'] < max_percent_mito, :];
        
    def scanpyQuickProcess(self):
        adata = self.data.copy();
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var['highly_variable']]
        sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')
        
        vr=list(adata.uns["pca"]["variance_ratio"]);
        vr = [round(x,5) for x in vr];
        dup=0;
        pre=""
        npcs=0;
        for i in range(len(vr)):
            if vr[i]==pre:
                dup+=1;
                if dup>=2:
                    npcs=i-2;
                    break;
            else:
                pre=vr[i];
                dup=0;

        npcs =40 if npcs>40 else npcs;
        npcs =4 if npcs<4 else npcs;
            
        print("n pcs: "+str(npcs));
        self.npcs=npcs;
        
        vr=list(adata.uns["pca"]["variance_ratio"]);
        vr = [round(x,4) for x in vr];
        dup=0;
        pre=""
        tsne_npcs=0;
        for i in range(len(vr)):
            if vr[i]==pre:
                dup+=1;
                if dup>=3:
                    tsne_npcs=i-1;
                    break;
            else:
                pre=vr[i];
                dup=0;

        tsne_npcs =30 if tsne_npcs>30 else tsne_npcs;
        tsne_npcs =4 if tsne_npcs<4 else tsne_npcs;
            
        print("tsne n pcs: "+str(tsne_npcs));
        self.tsne_npcs=tsne_npcs;
        
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=npcs);
        sc.settings.set_figure_params(dpi=80);
        sc.tl.louvain(adata)
        sc.tl.leiden(adata)
        sc.tl.umap(adata)
        sc.tl.phate(adata)
        sc.tl.tsne(adata,n_pcs=tsne_npcs);
        self.annoData=adata;
        print("done.")
        

        
    def showTsne(self,color='leiden'):
        sc.settings.set_figure_params(dpi=80)
        sc.pl.tsne(self.annoData, color=color)

    def showUmap(self,color='leiden'):
        sc.settings.set_figure_params(dpi=80)
        sc.pl.umap(self.annoData, color=color);
    
    def showPhate(self,color='leiden'):
        sc.settings.set_figure_params(dpi=80)
        sc.pl.phate(self.annoData, color=color);
    
    
    def creatNormalData(self):
        print("normalize....")
        matrix = self.data.to_df();
        diff= matrix.max(axis=1)-matrix.min(axis=1);
        matrix = matrix.sub(matrix.min(axis=1),axis=0).div(diff,axis=0).apply(lambda x: round(x.astype(float)*10000,2)  )
        self.ndata=matrix.T;
        
    
    
    def insertToDB(self,dbname='scDB',dbport=27017,dbhost="localhost",
                   adata="",mapType="tsne",
                   mapName='',
                   study="",
                   subjectid="",
                   disease="",
                   source="",
                   sample="",
                   comment="",
                   author="",
                  ):
        if mapName=="":
            print("mapName cannot empty");
            
        
        
        client = MongoClient(dbhost,dbport)
        db = client[dbname];
        
        coor=None;
        if adata =="":
            try:
                adata=self.annoData;
            except:
                print("annoData not exist.");
                return ""
        
            
        if mapType=='tsne':
            coor = adata.obsm["X_tsne"];
        elif mapType == "umap":
            coor = adata.obsm["X_umap"];
        elif mapType == "phate":
            coor = adata.obsm["X_phate"];
        else:
            print("no matched plot");
            return;
        
        cells =adata.obs.index
        metaDict=dict();
        leiden_val  = adata.obs["leiden"]
        for i in range(len(cells)):
            cellid = cells[i];
            metaDict[cells[i]]={"leiden":leiden_val[i]  };
        
        coordata=[];
        mapsamples=dict();
        samplesCluster=dict();
        
        for i in range(len(cells)):
            cell = cells[i];
            tx = str(coor[i][0]);
            ty= str(coor[i][1]);

            cell2=cell.split("-")[-1];
            if cell2.isnumeric():
                if cell2 in samplesCluster:
                    samplesCluster[cell2].append(i);
                else:
                    samplesCluster[cell2]=[i];
                    
            mapsamples[cell]=i;
            temp={"_id":cell,"x":round(float(tx),7),"y":round(float(ty),7),"order":i  };
            coordata.append(temp);
            
        if self.samples !="" and len(samplesCluster.keys())==len(self.samples):
            sampleClstrkeys=list(samplesCluster.keys());
            for i in sampleClstrkeys:
                i2=int(i)-1;
                newClstrName=self.samples[i2];
                samplesCluster[newClstrName]=samplesCluster.pop(i);
        
                
                    
        self.creatNormalData();
        
        exprhead=list(self.ndata.columns.values)
        
        exprhead_order=[];
        for i in exprhead:
            exprhead_order.append(mapsamples[i]);
            
        genes=list(self.ndata.index);
        print("start insert to db");
        
        mapinfo=dict();
        mapinfo["name"]= mapName;
        mapinfo["study"]=study;
        mapinfo["subjectid"]= subjectid ;
        mapinfo["tissue"]=sample;
        mapinfo["disease"]= disease;
        mapinfo["source"]= source ;
        mapinfo["comment"]= comment ;
        mapinfo["author"]= author ;
        

        #init end,start insert
        newmap = db.dataInfo.insert(mapinfo);
        newmap = str(newmap);
        exprCollection = "expr_"+newmap;
        for i in range(len(genes)):
            gene=genes[i].strip();
            expr=list(self.ndata.values[i]);
            newexpr=[0]*len(expr);
            for j in range(len(expr)):
                value=float(expr[j]);
                position=exprhead_order[j];
                if value !=0:
                    newexpr[position]=round(value,4);
            db[exprCollection].update({"_id":gene.upper().strip()},{"$set":{"normalize":newexpr}},upsert=True);
        
        metaCollection ='meta_'+newmap;
        for i in coordata:
            db[metaCollection].insert(i);
        
        colorlist = ["#EF4036","#907DBA" ,"#38B449","#F7931D","#F8ED31","#484B5A","#55B5E6","#F37E87","#70C38F" , "#3399FF","#0078AE", "#0000FF", "#A6E286", "#00AE3C",'#006400',"#F4C2C2","#FA6E79","#D1001C",'#660000','#FFD831','#FF8C00',"#AA6C39","#966FD6", "#B23CBF","#713E90", "#999999","#D1BEA8",'#54626F','#3B3C36']
        #metaDict[cells[i]]={"leiden":leiden_val[i]  };
        leidenclstr=dict();
        for i in metaDict:
            cell=i;
            clstrname=metaDict[i]["leiden"];
            if clstrname in leidenclstr:
                leidenclstr[clstrname].append(mapsamples[cell])
            else:
                leidenclstr[clstrname]=[mapsamples[cell]]
        
        indx=0;
        for i in leidenclstr:
            clstrType="leiden";
            tempcolor = colorlist[indx];
            db.cluster.insert({"mapid":ObjectId(newmap),"clstrType":clstrType,"clstrName":str(i),
                              "cells":leidenclstr[i],"color":tempcolor,
                               "x":"","y":"","label":False,"prerender":True
                              })
            indx +=1;
            if indx == len(colorlist):
                indx=0;
                
                
        indx=0;
        for i in samplesCluster:
            clstrType="samples";
            tempcolor = colorlist[indx];
            db.cluster.insert({"mapid":ObjectId(newmap),"clstrType":clstrType,"clstrName":str(i),
                              "cells":samplesCluster[i],"color":tempcolor,
                               "x":"","y":"","label":False,"prerender":True
                              })
            indx +=1;
            if indx == len(colorlist):
                indx=0;
                
        print("success")
        print("mapid: "+newmap);
        return newmap
