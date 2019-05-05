import os, sys, csv,json,datetime,time,math,scipy.stats, collections;
import re

import numpy as np;
import scipy as sp;
import math;
import linecache

from scipy.spatial import distance 
from scipy.stats import ranksums;

import pymongo;
from bson.objectid import ObjectId
from pymongo import MongoClient


from bson.json_util import dumps

import sklearn;


#import copy
import random
#import operator;

import multiprocessing as mp


client = MongoClient('localhost');

databaseName="scDB"



db = client[databaseName];



def loginVerify(username,password):

	user = db.userInfo.find({"_id":username,"password":password})
	usercount = user.count();

	userrole ="";
	if usercount ==1:
		userrole = user[0]["role"];
	else:
		userrole ="";

	data = dict();
	data["role"]=userrole;
	data["count"]=usercount;

	return data;



def checkDataLength(mapid):

	return "";



def getMapInfoBySampleId(sampleid):
	sampleid = ObjectId(sampleid);
	data = db.dataInfo.find_one( {"_id":sampleid},{"name":1})
	info =dict();
	info["mapname"]=data["name"];
	
	return info;

def getMapDataBySampleId(sampleid):
	
	tsneMetaCollection = "meta_"+sampleid;

	tsneMap = db[tsneMetaCollection].find({}).sort([("order", pymongo.ASCENDING)]);

	tsneMap = list(tsneMap);

	return tsneMap;

def getAllSampleInfo(userid):
	data = db.dataInfo.find({},{"_id":1,"name":1,"study":1,"subjectid":1,"tissue":1,"disease":1,"source":1,"comment":1});
	resdata=[];
	for i in data:
		i["_id"]=str(i["_id"]);
		resdata.append(i);
	return resdata;



def listClusters(mapid):
	mapid = ObjectId(mapid);
	clstrs = db.cluster.distinct("clstrType",{"mapid":mapid});
	clstrs = list(clstrs);
	return clstrs;


def getClusterInfo2(mapid,clstrType):

	mapid = ObjectId(mapid);
	clstrs = db.cluster.find({"mapid":mapid,"clstrType":clstrType},{"_id":1,"clstrName":1,"color" :1, "x":1, "y" : 1, "label" : 1, "prerender" : 1,"marks":1,"negmarks":1,"cells":1});

	res=[];
	for i in clstrs:
		i["_id"] = str(i["_id"]);
		res.append(i)

	return res;

def getClusterInfo3(mapid,clstrType):

	mapid = ObjectId(mapid);
	clstrs = db.cluster.find({"mapid":mapid,"clstrType":clstrType},{"_id":1,"clstrName":1});

	res=[];
	for i in clstrs:
		i["_id"] = str(i["_id"]);
		res.append(i)

	return res;

def getClusterInfo(sampleid):


	sampleid = ObjectId(sampleid);
	clstrs = db.cluster.find({"mapid":sampleid},{"_id":1,"clstrName":1,"clstrType":1,"color" :1, "x":1, "y" : 1, "label" : 1, "prerender" : 1,"marks":1,"negmarks":1});

	resclstrs=dict();
	for i in clstrs:
		clstrtype = i["clstrType"];
		idstr=str(i["_id"]);
		i["_id"]=idstr;
		if clstrtype in resclstrs:
			resclstrs[clstrtype].append(i);
		else:
			resclstrs[clstrtype]=[i];

	return resclstrs;


def getExprdataByGene(sampleid,gene):
	
	genexpr = db["expr_"+sampleid].find_one({"_id":gene},{"normalize":1});

	if genexpr == None:
		return None;
	countexpr= genexpr["normalize"];

	res =dict();

	for i in range(len(countexpr)):
		if countexpr[i] >0:
			res[i]=countexpr[i];


	return res;


def queryClstrType():
	clstrTypes = db.clusterType.distinct("_id");
	clstrTypes = list(clstrTypes);
	return clstrTypes;




def deleteMapById(mapid):
	db.dataInfo.remove({"_id":ObjectId(mapid)});
	db.cluster.remove({"mapid":ObjectId(mapid)});
	db.drop_collection("expr_"+mapid)
	db.drop_collection("meta_"+mapid)


	return "success"



def updateMap(mapid,name,subjectid,source,study,disease,tissue):
	mapid = ObjectId(mapid);
	res = db.dataInfo.update({"_id":mapid},{"$set":{ "name":name,"subjectid":subjectid,"source": source,"study":study,"disease":disease ,"tissue":tissue   }},upsert=False);
	#['n', 'nModified', 'ok', 'updatedExisting']

	return "success"



def getMapInfoByDiseaseCategory(diseaseCategory):

	if diseaseCategory == "ALL":
		res = db.dataInfo.find();
		
	"""
	elif diseaseCategory == "CNS":
		res = db.dataInfo.find({"disease":"CNS"});


	elif diseaseCategory == "Immunology":
		
		res = db.dataInfo.find({ "$or":[ 
			{"disease":{"$regex":"icd","$options":"$i"}},
			
			{"disease":{"$regex":"uc","$options":"$i"}},
		]  });

		

	elif diseaseCategory == "MetabolicDisease":
		res = db.dataInfo.find({ "$or":[ {"study": {"$regex":"^Chromium_20181130","$options":"$i"}},{"study": "CircuitObesity"}]  });
		

	elif diseaseCategory == "Respiratory":
		res = db.dataInfo.find({"disease":"COPD"});

	elif diseaseCategory == "CancerImmunology":
		res = db.dataInfo.find({ "$or":[ 
			{"study": {"$regex":"cancer","$options":"$i"}},
			{"subjectid":  {"$regex":"cancer","$options":"$i"} },
			{"disease":  {"$regex":"cancer","$options":"$i"} },
			{"tissue":  {"$regex":"cancer","$options":"$i"} },
			{"source":  {"$regex":"cancer","$options":"$i"} },
			{"disease":{"$regex":"Melanoma","$options":"$i"}},
			{"disease":{"$regex":"HCC","$options":"$i"}},
			{"disease":{"$regex":"NSCLC","$options":"$i"}},
			{"disease":{"$regex":"Lung Adeno Ca","$options":"$i"}},
		]  });

	elif diseaseCategory == "CKD":
		res = db.dataInfo.find({"disease":"CKD"});
	elif diseaseCategory == "CellAtlas":
		res = db.dataInfo.find({ "$or":[ {"study": {"$regex":"pbmc","$options":"$i"}},
			{"subjectid":  {"$regex":"pbmc","$options":"$i"} },
			{"tissue":  {"$regex":"pbmc","$options":"$i"} },
			{"tissue":  {"$regex":"E-MTAB-6678","$options":"$i"} },
			{"source":  {"$regex":"Public Regeneron Zhao Q","$options":"$i"} },
			{"disease":  {"$regex":"Healthy","$options":"$i"} },
			{"disease":  {"$regex":"Normal","$options":"$i"} },
		]  });
	"""
	res2 = dict();
	for i in res:
		key=  i["disease"]+i["study"]+i["source"]+ i["subjectid"];
		
		if key in res2:
			res2[key].append({"_id":str(i["_id"]),"disease":i["disease"],"study":i["study"],"source":i["source"],"tissue":i["tissue"],"subjectid":i["subjectid"],"name":i["name"] });
		else:
			res2[key]=[{"_id":str(i["_id"]),"disease":i["disease"],"study":i["study"],"source":i["source"],"tissue":i["tissue"],"subjectid":i["subjectid"],"name":i["name"] }]


	res3=[];
	for i in res2:
		datas= res2[i];
		data=datas[0];
		common={"disease":data["disease"],"study":data["study"],"source":data["source"] ,"subjectid":data["subjectid"]};
		temp=[];
		for j in datas:
			temp.append({"cid":j["_id"],"sample":j["tissue"],"name":j["name"]});

		common["samples"]=temp
		res3.append(common)

	return res3;



def getGeneSearchPlotDataBycellType(gene,mapid,clstrType):

	genexpr = db["expr_"+mapid].find_one({"_id":gene},{"normalize":1});
	if genexpr == None:
		return None;
	countexpr= genexpr["normalize"];
	clusters = db.cluster.find({"mapid":ObjectId(mapid),"clstrType":clstrType},{"cells":1,"clstrName" : 1, "_id":1,"color":1});
	res=[]
	for i in clusters:
		cid = str(i["_id"]);
		cells=i["cells"]
		clstrlen = len(cells);
		nonzeros=[];
		for pos in cells:
			expr_val = countexpr[pos];
			if expr_val >0:
				nonzeros.append(expr_val);
		if len(nonzeros)==0:
			nonzeros_mean = 0;
			nonzeros_median = 0;
			nonzeros_1percentile=0;
			nonzeros_3percentile=0;
			nonzeros_perc = 0;
			nonzeros_min = 0;
			nonzeros_max = 0;
			p100=0;
			p0=0;
		else:
			nonzeros_mean = np.mean(nonzeros);
			
			nonzeros_1percentile = np.percentile(nonzeros,25);
			nonzeros_median = np.percentile(nonzeros,50);
			nonzeros_3percentile = np.percentile(nonzeros,75);
			nonzeros_perc = len(nonzeros)/clstrlen;
			#print(nonzeros_1percentile);
			#print(nonzeros_median)
			#print(nonzeros_3percentile);
			
			iqr = nonzeros_3percentile-nonzeros_1percentile;
			nonzeros_min = nonzeros_1percentile-1.5*iqr;
			nonzeros_max = nonzeros_3percentile+1.5*iqr;

			p100=np.percentile(nonzeros,100);
			p0=np.percentile(nonzeros,0);
			if nonzeros_min < p0:
				nonzeros_min=p0

			if nonzeros_max > p100:
				nonzeros_max=p100;

		res.append({"name":i["clstrName"],"color":i["color"],"cid":cid,"mean":nonzeros_mean,"median":nonzeros_median,"perc":nonzeros_perc,"1q":nonzeros_1percentile,"3q":nonzeros_3percentile,"min":nonzeros_min,"max":nonzeros_max,"p100":p100,"p0":p0,"cells":cells})

	return res;

def getGeneSearchPlotData(gene,sampleid):
	
	genexpr = db["expr_"+sampleid].find_one({"_id":gene},{"normalize":1});

	if genexpr == None:
		return None;
	countexpr= genexpr["normalize"];

	clusters = db.cluster.find({"mapid":ObjectId(sampleid)},{"cells":1,"clstrName" : 1, "clstrType":1,"_id":1,"color":1});

	res=[]
	for i in clusters:
		
		cid = str(i["_id"]);
		cells=i["cells"]
		clstrlen = len(cells);
		nonzeros=[];
		for pos in cells:
			expr_val = countexpr[pos];
			if expr_val >0:
				nonzeros.append(expr_val);

		if len(nonzeros)==0:
			nonzeros_mean = 0;
			nonzeros_median = 0;
			nonzeros_1percentile=0;
			nonzeros_3percentile=0;
			nonzeros_perc = 0;
			nonzeros_min = 0;
			nonzeros_max = 0;
			p100=0;
			p0=0;
		else:
			nonzeros_mean = np.mean(nonzeros);
			
			nonzeros_1percentile = np.percentile(nonzeros,25);
			nonzeros_median = np.percentile(nonzeros,50);
			nonzeros_3percentile = np.percentile(nonzeros,75);
			nonzeros_perc = len(nonzeros)/clstrlen;
			#print(nonzeros_1percentile);
			#print(nonzeros_median)
			#print(nonzeros_3percentile);
			
			iqr = nonzeros_3percentile-nonzeros_1percentile;
			nonzeros_min = nonzeros_1percentile-1.5*iqr;
			nonzeros_max = nonzeros_3percentile+1.5*iqr;

			p100=np.percentile(nonzeros,100);
			p0=np.percentile(nonzeros,0);
			if nonzeros_min < p0:
				nonzeros_min=p0

			if nonzeros_max > p100:
				nonzeros_max=p100


		res.append({"name":i["clstrName"],"color":i["color"],"ctype":i["clstrType"],"cid":cid,"mean":nonzeros_mean,"median":nonzeros_median,"perc":nonzeros_perc,"1q":nonzeros_1percentile,"3q":nonzeros_3percentile,"min":nonzeros_min,"max":nonzeros_max,"p100":p100,"p0":p0})
		#break;

	return res;





def getExprNormailizedataByGene(sampleid,gene):
	
	genexpr = db["expr_"+sampleid].find_one({"_id":gene},{"normalize":1});

	if genexpr == None:
		return None;
	countexpr= genexpr["normalize"];

	res =[];
	for i in range(len(countexpr)):
		if countexpr[i] >0:
			res.append(i);
	
	return res;




def listExistsGenes(sampleid,genes):

	fitGenes= db["expr_"+sampleid].distinct("_id",{"_id":{"$in":genes}})

	return fitGenes;

def listExistsGenesRegex(sampleid,geneRegex):

	fitGenes= db["expr_"+sampleid].distinct("_id",{"_id":{"$regex":"^"+geneRegex,"$options": "i" }})
	
	return fitGenes;


def savecluster(sampleid,name,ctype,cells,comment,marks,negmarks):

	sampleid = ObjectId(sampleid);

	color =getRandomColor();

	clstrcount = db.cluster.find({"mapid":sampleid,"clstrName":name,"clstrType":ctype}).count();

	if clstrcount == 0:
		comment = comment.strip();
		db.cluster.insert_one({"mapid":sampleid,"clstrName":name,"clstrType":ctype,"color":color,"cells":cells,"comment":comment,"x" : "", "y" : "", "label" : False, "prerender" : True,"marks":marks,"negmarks":negmarks})
		clstr = db.cluster.find_one({"mapid":sampleid,"clstrName":name,"clstrType":ctype});

		clstr_id=str(clstr["_id"])
		return {"status":"success","cid":clstr_id,"color":color}

	else:

		return {"status":"failed"}



def getClusterCellsById(clstrid):
	clstrid = ObjectId(clstrid);
	clstrCells = db.cluster.find_one({"_id":clstrid},{"cells":1});
	clstrCells= clstrCells["cells"];
	return clstrCells;


def contrastCellsVsClstr(sampleid,cells,clstr,contrastId):
	cell2 = getClusterCellsById(clstr);
	cell1= np.array(cells,dtype='i');
	data = contrast(sampleid,cell1,cell2);
	db2=MongoClient('localhost')[databaseName];
	db2.contrastResult.update_one({"_id":contrastId},{"$set":{"p":data["p"],"n":data["n"],"done":True} })

	return "";


def queryClstrCellsAndLabelByCid(cid):
	cid = ObjectId(cid);
	clstr = db.cluster.find_one({"_id":cid},{"cells":1,"label":1,"x":1,"y":1,"clstrName":1});
	
	res=dict();
	res["cellids"]=clstr["cells"];
	res["name"] = clstr["clstrName"];
	res["labeled"]=clstr["label"];
	res["x"] = clstr["x"];
	res["y"]= clstr["y"];

	return res;



def getClusterClassification(clstrtype):

	data = db[clstrtype].find({});

	return list(data);


def getExprPosCountsByGene(sampleid,g1,g2):
	genexpr1 = db["expr_"+sampleid].find_one({"_id":g1},{"normalize":1})["normalize"];
	genexpr2 = db["expr_"+sampleid].find_one({"_id":g2},{"normalize":1})["normalize"];

	c1=0;
	c2=0;
	c3=0;
	d1=[];
	d2=[];
	d3=[];
	for i in range(len(genexpr1)):
		v1 = genexpr1[i];
		v2 = genexpr2[i];

		if v1>0 and v2 >0:

			c3+=1;
			d3.append(i);
		else:
			if v1>0:
				c1+=1;
				d1.append(i);

			if v2>0:
				c2+=1;
				d2.append(i);

	return {"g":[g1,g2,"intersc"],"c":[c1,c2,c3],"d1":d1,"d2":d2,"d3":d3}



def getExprPosCountByGene(sampleid,gene):
	genexpr = db["expr_"+sampleid].find_one({"_id":gene},{"normalize":1});

	if genexpr == None:
		return None;
	countexpr= genexpr["normalize"];

	count=0;
	for i in countexpr:
		if i >0:
			count+=1;


	return count;


def updateClusterColor(clstrid,color):
	res = db.cluster.update_one({"_id":ObjectId(clstrid)},{"$set":{"color": color}});

	return "success";

def updateClusterPostition(clstrid,x,y):
	res = db.cluster.update_one({"_id":ObjectId(clstrid)},{"$set":{"x":float(x),"y":float(y),"label":True}});
	return "success";

def updateClusterName(clstrid,name):
	res = db.cluster.update_one({"_id":ObjectId(clstrid)},{"$set":{"clstrName":name}});
	return "success";

def deleteCluster(clstrid):
	res = db.cluster.remove({"_id":ObjectId(clstrid)});

	return "success";


def updateClusterMarks(clstrid,marks):
	res = db.cluster.update_one({"_id":ObjectId(clstrid)},{"$set":{"marks":marks}});

	return "success";


def updateClusterNegMarks(clstrid,marks):
	res = db.cluster.update_one({"_id":ObjectId(clstrid)},{"$set":{"negmarks":marks}});

	return "success";


def updateClusterIsPreRender(clstrid,val):
	if val =='T':
		val =True;
	elif val =="F":
		val = False;
	res = db.cluster.update_one({"_id":ObjectId(clstrid)},{"$set":{"prerender":val}});
	return "success";





def contrast(sampleid,cells1,cells2):

	xlen = len(cells1);
	ylen = len(cells2);

	p=dict();
	n=dict();
	allexpr = db["expr_"+sampleid].find({},{"_id":1,"normalize":1});
	for i in allexpr:
		g = i["_id"];
		expr = i["normalize"];
		x=[];
		y=[];

		xpos=0;
		ypos=0;

		for j in cells1:
			x.append(expr[j])
			if expr[j]>0:
				xpos+=1;


		for j in cells2:
			y.append(expr[j]);
			if expr[j]>0:
				ypos+=1;


		percx = xpos/xlen;
		percy = ypos/ylen;

		
		if percx> 0.95 and percy > 0.95:
			continue;
		

		ranksumsres = scipy.stats.ranksums(x,y);	
		statics = ranksumsres[0];
		pval = ranksumsres[1];
			
		if pval < 0.01:
			if statics >=0:
				p[g]=pval;
			else:
				n[g]=pval;

	p = sorted(p.items(), key=lambda kv: kv[1] );
	n =	sorted(n.items(), key=lambda kv: kv[1] );

	p2=[];
	n2=[];
	for i in p:
		p2.append(i[0]);
	p=None;
	for i in n:
		n2.append(i[0]);
	n=None;
	return {"p":p2,"n":n2}




def getContrastResult(resultid):

	resultid = ObjectId(resultid);

	data = db.contrastResult.find_one({"_id":resultid});
	
	done = data["done"];
	if done:
		p=data["p"];
		n=data["n"];
		db.contrastResult.remove({"_id":resultid});
		return {"done":done,"p":p,"n":n}
	else:
		return {"done":done}






def doContrast(sampleid,cells,clstr,contrastModel):
	currentTime = time.time()
	contrastId = db.contrastResult.insert({"startTime":currentTime,"done":False}); 
	returncontrastId = str(contrastId);


	if contrastModel == "contrastwithrest":
		p1 = mp.Process(target=contrastwithrest,args=(sampleid,cells,contrastId) )
		p1.start();

	if contrastModel == "contrastCellsVsClstr":

		p1 = mp.Process(target=contrastCellsVsClstr,args=(sampleid,cells,clstr,contrastId) )
		p1.start();


	return returncontrastId;


def contrastwithrest(sampleid,cells,contrastId):	
	#db.contrastResult.insert_one()
	db2=MongoClient('localhost')[databaseName];

	#cells=np.array(cells,dtype="i");
	cellsdict= dict();
	xlen=0
	for i in cells:
		xlen+=1;
		cellsdict[int(i)]=None;



	p=dict();
	n=dict();

	allexpr = db["expr_"+sampleid].find({},{"_id":1,"normalize":1});
	for i in allexpr:
		g = i["_id"];
		expr = i["normalize"];
		#filter

		totallen = len(expr)
		ylen = totallen-xlen;



		xpos=0;
		ypos=0;
		x=[];
		y=[];
		for j in range(len(expr)):
			if j in cellsdict:
				x.append(expr[j]);
				if expr[j]>0:
					xpos+=1;
			else:
				y.append(expr[j]);
				if expr[j]>0:
					ypos+=1;



		percx = xpos/xlen;
		percy = ypos/ylen;

		
		if percx> 0.95 and percy > 0.95:
			continue;


		ranksumsres = scipy.stats.ranksums(x,y);
		
		statics = ranksumsres[0];
		pval = ranksumsres[1];

		
		
		if pval < 0.01:
			if statics >=0:
				p[g]=pval;
			else:
				n[g]=pval;

	p = sorted(p.items(), key=lambda kv: kv[1] );
	n =	sorted(n.items(), key=lambda kv: kv[1] );

	p2=[];
	n2=[];
	for i in p:
		p2.append(i[0]);
	p=None;
	for i in n:
		n2.append(i[0]);
	n=None;
	
	db2.contrastResult.update_one({"_id":contrastId},{"$set":{"p":p2,"n":n2,"done":True}} );

	#return {"p":p2,"n":n2}
	return "";

def runRanksums(sampleid,arr,compareTargets):

	if compareTargets =='all':
		pass;
	else:
		print(compareTargets);





	return ""


def contrastGeneSearch(gene,cells1,cells2,sampleid,name1,name2):


	genexpr = db["expr_"+sampleid].find_one({"_id":gene},{"normalize":1});

	if genexpr == None:
		return None;
	countexpr= genexpr["normalize"];

	maxval =0;

	exprdict =dict();
	for i in range(len(countexpr)):
		if countexpr[i] >0:
			exprdict[i]=countexpr[i];
			if countexpr[i] > maxval:
				maxval	= countexpr[i];
	
	res=[];

	for i in [1,2]:
		if i ==1:
			cells=cells1;
			name=name1;
			color="orange";
		elif i ==2:
			cells=cells2;
			name=name2;
			color="steelblue";

		clstrlen = len(cells);
		nonzeros=[];
		for pos in cells:
			expr_val = countexpr[pos];
			if expr_val >0:
				nonzeros.append(expr_val);

		if len(nonzeros)==0:
			nonzeros_mean = 0;
			nonzeros_median = 0;
			nonzeros_1percentile=0;
			nonzeros_3percentile=0;
			nonzeros_perc = 0;
			nonzeros_min = 0;
			nonzeros_max = 0;
			p100=0;
			p0=0;
		else:
			nonzeros_mean = np.mean(nonzeros);
			
			nonzeros_1percentile = np.percentile(nonzeros,25);
			nonzeros_median = np.percentile(nonzeros,50);
			nonzeros_3percentile = np.percentile(nonzeros,75);
			nonzeros_perc = len(nonzeros)/clstrlen;
			#print(nonzeros_1percentile);
			#print(nonzeros_median)
			#print(nonzeros_3percentile);
			
			iqr = nonzeros_3percentile-nonzeros_1percentile;
			nonzeros_min = nonzeros_1percentile-1.5*iqr;
			nonzeros_max = nonzeros_3percentile+1.5*iqr;

			p100=np.percentile(nonzeros,100);
			p0=np.percentile(nonzeros,0);
			if nonzeros_min < p0:
				nonzeros_min=p0

			if nonzeros_max > p100:
				nonzeros_max=p100


		res.append({"name":name,"color": color,"mean":nonzeros_mean,"median":nonzeros_median,"perc":nonzeros_perc,"1q":nonzeros_1percentile,"3q":nonzeros_3percentile,"min":nonzeros_min,"max":nonzeros_max,"p100":p100,"p0":p0})
	

	return res,exprdict,maxval;



#import multiprocessing as mp;
#p1=mp.Process(target=runWilcoxon,args=(arr))
#p1.start();
#p1.join();


def getClusterRestCells(sampleid,cells):
	
	cellorders = db["meta_"+sampleid].find({},{"order":1});
	
	cellsdict=dict();
	for i in cells:
		cellsdict[i]=None;


	res=[];
	for i in cellorders:
		if i["order"] not in cellsdict:
				res.append(i["order"]);

	return res;




def strarrayToIntarray(cellstr):
	cells = np.array(cellstr.split(","),dtype='i');
	return list(cells);



def getRandomColor():
	r = lambda: random.randint(0,255);
	color = '#%02X%02X%02X' % (r(),r(),r());

	return color;






def getClusterNameById(cid):
	clstrid = ObjectId(cid);
	clstrName = db.cluster.find_one({"_id":clstrid},{"clstrName":1});
	clstrName= clstrName["clstrName"];
	return clstrName;



def getAllClusterStudies():

	studies = db.dataInfo.aggregate([
		{"$group":{
			 "_id":"$study" ,
			 "tissues":{"$addToSet":"$tissue"} 
		}}
	]);


	studies =list(studies);

	return studies;



def getAllTissueByStudies(study):
	data = db.dataInfo.find({"study":study});




def getAllClusterTypesByStudyAndTissues(study,tissues):


	mapids = db.dataInfo.distinct("_id",{"study":study,"tissue":{"$in":tissues}});
	data = db.cluster.aggregate([
		{"$match":{"mapid":{"$in" :mapids}}},

		{"$project":{"clstrName":1,"clstrType":1}},
		{"$group":{"_id":"$clstrType","clstrName":{"$addToSet":"$clstrName"}}},
		{"$project":{"clstrType":"$_id","clstrName":1,"_id":0}},


	])

	datadict=dict();
	for i in data:
		datadict[i["clstrType"]]=i["clstrName"];


	return datadict;



def queryComparePlotData(study,tissues,gene,cluster,clusterType):
	
	gene = gene.upper();

	maps = db.dataInfo.find({"study":study,"tissue":{"$in":tissues}},{"_id":1,"tissue":1,"subjectid":1,"name":1});

	mapids=[];

	result = dict();

	for i in maps:
		mapids.append(i["_id"]);
		result[str(i["_id"])]={"tissue":i["tissue"],"subjid":i["subjectid"],"name":i["name"]};

	clusterinfo = db.cluster.find( {"mapid":{"$in" :mapids},"clstrName":cluster,"clstrType":clusterType },{"mapid":1,"cells":1}   )
	for i in clusterinfo:

		mapid = str(i["mapid"]);
		cells = i["cells"];
		expr = db["expr_"+mapid].find_one({"_id":gene},{"gene":1,"normalize":1});
		if expr is not None:
			normVal = expr["normalize"];

			nonzeros=[];
			allexpr =[]
			for pos in cells:
				expr_val = normVal[pos];
				if expr_val >0:
					nonzeros.append(expr_val);

				allexpr.append(expr_val);

			if len(nonzeros)==0:
				nonzeros_mean = 0;
				nonzeros_median = 0;
				nonzeros_1percentile=0;
				nonzeros_3percentile=0;
				nonzeros_perc = 0;
				nonzeros_min = 0;
				nonzeros_max = 0;
				p100=0;
				p0=0;
			else:
				nonzeros_mean = np.mean(nonzeros);
				
				nonzeros_1percentile = np.percentile(nonzeros,25);
				nonzeros_median = np.percentile(nonzeros,50);
				nonzeros_3percentile = np.percentile(nonzeros,75);
				nonzeros_perc = len(nonzeros)/len(allexpr);
				
				iqr = nonzeros_3percentile-nonzeros_1percentile;
				nonzeros_min = nonzeros_1percentile-1.5*iqr;
				nonzeros_max = nonzeros_3percentile+1.5*iqr;

				p100=np.percentile(nonzeros,100);
				p0=np.percentile(nonzeros,0);
				if nonzeros_min < p0:
					nonzeros_min=p0

				if nonzeros_max > p100:
					nonzeros_max=p100


			result[mapid]["expr"]=allexpr;
			result[mapid]["mean"]=nonzeros_mean;
			result[mapid]["median"]=nonzeros_median;
			result[mapid]["perc"]=nonzeros_perc;
			result[mapid]["1q"]=nonzeros_1percentile;
			result[mapid]["3q"]=nonzeros_3percentile;
			result[mapid]["min"]=nonzeros_min;
			result[mapid]["max"]=nonzeros_max;
			result[mapid]["p100"]=p100;
			result[mapid]["p0"]=p0;

	result2=dict();

	for i in result:
		subj = result[i]["subjid"];
		if "expr" in result[i]:
			if subj in result2:
				result2[subj].append(result[i]);
				
			else:
				result2[subj]=[result[i]];

	return result2;

def getBarcodes(cellids,mapid):
	cellids2 =[];
	for i in cellids:
		cellids2.append(int(i));
	barcodes = db["meta_"+mapid].distinct("_id",{"order":{"$in":cellids2} } );

	return barcodes;


def getCellsByStudyAndTissueAndclstrtypeAndClstr(study,tissue,clusterType,cluster):

	study = str(study).strip();
	clusterType = str(clusterType).strip();
	cluster = str(cluster).strip();
	tissue = str(tissue).strip();

	maps = db.dataInfo.find({"study":study,"tissue": tissue},{"_id":1,"subjectid":1,"tissue":1});
	
	result=dict();
	mapids=[]
	for i in maps:
		result[str(i["_id"])]={"name": str(i["subjectid"])+"__"+str(i["tissue"])};
		mapids.append(i["_id"]);
	clstrs = db.cluster.aggregate([
		{"$match":{"mapid":{"$in":mapids},"clstrType":clusterType,"clstrName":cluster}},
		{"$project":{"_id":0,"mapid":1,"cells":1}},
	]);
	for i in clstrs:
		cells = db["meta_"+str(i["mapid"])].find({"order":{"$in":i["cells"]} },{"_id":1,"order":1});

		result[str(i["mapid"])]["cells"]=list(cells);

	"""
	result2=dict();
	for i in result:
		mapid = i;
		name = result[i]["name"].replace(" ","_");
		cells = result[i]["cells"];

		for c in cells:
			cname = name+"_"+c["_id"];
			order = c["order"];

			result2[cname]={"id":mapid,"odr":order};

	return result2;

	"""
	return result;

#def getALLCellTypeByStudy(study):






"""


def getPatientProfilerDataByGene(study,ppname,gene,isagg):

	pp=patientProfilerDB.meta.find_one({"scope":study,"name":ppname});

	tumorSamples = exprMatrixDB[study+"_meta"].distinct("order",{"grouplist.category":"TUMOR OR NORMAL","grouplist.group":"TUMOR"})
	normalSamples = exprMatrixDB[study+"_meta"].distinct("order",{"grouplist.category":"TUMOR OR NORMAL","grouplist.group":"NORMAL"})
	
	if len(tumorSamples) == 0 :
		tumorSamples = exprMatrixDB[study+"_meta"].distinct("order",{"grouplist.category":"TUMOR.OR.NORMAL","grouplist.group":"TUMOR"})
		normalSamples = exprMatrixDB[study+"_meta"].distinct("order",{"grouplist.category":"TUMOR.OR.NORMAL","grouplist.group":"NORMAL"})
	


	test_samples = pp["test_samples"];
	control_samples = pp["control_samples"];

	test_study = pp["test_study"];
	control_study = pp["control_study"];

	test_samples = exprMatrixDB[test_study+"_meta"].distinct("order",{"_id":{"$in":test_samples }});

	control_samples = exprMatrixDB[control_study+"_meta"].distinct("order",{"_id":{"$in": control_samples }});

	exprstr="$expr"
	if "isRNAseq":
		exprstr='$log2expr';

	aggEleAtarr=[];
	for i in tumorSamples:
		aggEleAtarr.append({"$arrayElemAt":[exprstr,i]});

	tumordt = exprMatrixDB[study].aggregate([
		{"$match":{"_id":gene}},
		{"$project":{"_id":0,"expr": aggEleAtarr}},
		{"$unwind":"$expr"},
		{"$sort":{"expr":1}}
	])

	tumor=[];
	for i in tumordt:
		tumor.append(round( i["expr"],3))


	aggEleAtarr=[];
	for i in normalSamples:
		aggEleAtarr.append({"$arrayElemAt":[exprstr,i]});

	normaldt = exprMatrixDB[study].aggregate([
		{"$match":{"_id":gene}},
		{"$project":{"_id":0,"expr": aggEleAtarr}},
		{"$unwind":"$expr"},
		{"$sort":{"expr":1}}

	])
	normal=[];
	for i in normaldt:
		normal.append(round( i["expr"],3))

	aggEleAtarr=[];
	for i in test_samples:
		aggEleAtarr.append({"$arrayElemAt":[exprstr,i]});

	testdt = exprMatrixDB[test_study].aggregate([
		{"$match":{"_id":gene}},
		{"$project":{"_id":0,"expr": aggEleAtarr}},
		{"$unwind":"$expr"},
		{"$sort":{"expr":1}}

	]);
	test = [];
	for i in testdt:
		test.append(round( i["expr"],3))

	aggEleAtarr=[];
	for i in control_samples:
		aggEleAtarr.append({"$arrayElemAt":[exprstr,i]});

	ctrldt = exprMatrixDB[control_study].aggregate([
		{"$match":{"_id":gene}},
		{"$project":{"_id":0,"expr": aggEleAtarr}},
		{"$unwind":"$expr"},
		{"$sort":{"expr":1}}
	]);

	ctrl = [];
	for i in ctrldt:
		ctrl.append(round( i["expr"],3));

	return {"s1":{"name":"Tumor" ,"values":tumor},
			"s2":{"name":"Normal","values":normal},
			"d":{"name":test_study ,"values":test},
			"n":{"name":control_study ,"values":ctrl}}


"""
