from django.shortcuts import render
import json
import os;
from django.http import HttpResponse, JsonResponse,HttpResponseRedirect

# Create your views here.

from . import service;




# C""reate your views here.
def index(request):

	return render(request,"index/index.html");


def neuralNetwork(request):

	return render(request,"index/neuralNetwork.html");

def rendermap(request,mapid):

	return render(request,"index/scmap.html",{"mapid":mapid} );


def compareplot(request):

	return render(request,"index/compareplot.html" );


def administration(request):

	return render(request,"index/admin.html");


def renderdatatable(request,diseaseCategory):

	return render(request,"index/mapTable.html",{"diseaseCategory":diseaseCategory})

def getMapInfoByDiseaseCategory(request):

	diseaseCategory = request.POST.get("diseaseCategory").strip();

	data = service.getMapInfoByDiseaseCategory(diseaseCategory);

	return JsonResponse({"data":data });

def querybarcodes(request):

	cellids = request.POST.get("cellids");
	cellids = cellids.split(",");
	sampleid = request.POST.get("sampleid");
	
	data = service.getBarcodes(cellids,sampleid);
	return JsonResponse({"data":data });

def updateMap(request):
	mapid = request.POST.get("_id").strip();
	name = request.POST.get("name").strip();
	subjectid = request.POST.get("subjectid").strip();
	source = request.POST.get("source").strip();
	study = request.POST.get("study").strip();
	disease = request.POST.get("disease").strip();
	tissue = request.POST.get("tissue").strip();
	res =  service.updateMap(mapid,name,subjectid,source,study,disease,tissue);
	

	return JsonResponse({"status":res });

def userLogin(request):
	username = request.POST.get("username");
	password = request.POST.get("password");
	data = service.loginVerify(username,password);
	
	status = "";
	if data["count"] ==1:
		status="success";
	else:
		status ="failed"

	return JsonResponse({"status":status,"role":data["role"]});


def getClstrsByTypeAndMapid2(request):
	mapid = request.POST.get("mapid");
	clstrType = request.POST.get("clstrType");
	clstrs = service.getClusterInfo3(mapid,clstrType);

	return JsonResponse({"res":clstrs});

def getClstrsByTypeAndMapid(request):
	mapid = request.POST.get("mapid");
	clstrType = request.POST.get("clstrType");
	clstrs = service.getClusterInfo2(mapid,clstrType);

	return JsonResponse({"res":clstrs});

def getMapDataByMapId(request):
	mapid = request.POST.get("mapid");

	data = service.getMapDataBySampleId(mapid);
	clstrInfo = service.listClusters(mapid);

	mapinfo = service.getMapInfoBySampleId(mapid);

	return JsonResponse({"tsneData":data,"clstr":clstrInfo,"info":mapinfo});

def getMapDataBySampleId(request):
	
	sampleid = request.POST.get("sampleid");

	data = service.getMapDataBySampleId(sampleid);
	clstrInfo = service.getClusterInfo(sampleid);

	mapinfo = service.getMapInfoBySampleId(sampleid);

	return JsonResponse({"tsneData":data,"clstr":clstrInfo,"info":mapinfo});

def getGeneSearchPlotData(request):
	gene = request.POST.get("gene");
	spid = request.POST.get("spid");
	data = service.getGeneSearchPlotData(gene,spid);
	
	return JsonResponse({"data": data});


def queryClstrType(request):
	data = service.queryClstrType();
	return JsonResponse({"data": data}); 


def deleteMapById(request):
	
	mapid = request.POST.get("mapid");
	service.deleteMapById(mapid);

	return JsonResponse({"success": "success"});


def getGeneSearchPlotDataByclstrType(request):

	mapid = request.POST.get("mapid");
	gene = request.POST.get("gene");
	clstrType = request.POST.get("clstrType");

	data = service.getGeneSearchPlotDataBycellType(gene,mapid,clstrType);

	return JsonResponse({"data":data});


def genelistSearch(request):
	sampleid = request.POST.get("sampleid");
	genestr = request.POST.get("genestr");

	gene = genestr.replace(" ","");
	gene = gene.split(",");
	
	gene2=[];
	for i in gene:
		if len(i) > 0:
			gene2.append(i.upper());


	data="";
	geneCount=0;
	gene="";
	if len(gene2) >1:
		data = service.listExistsGenes(sampleid,gene2);
		geneCount = len(data);
		if geneCount == 1:
			gene = data[0];
			data = service.getExprdataByGene(sampleid,data[0]);

		elif geneCount ==2:
			g1 = data[0];
			g2 = data[1];

			data = service.getExprPosCountsByGene(sampleid,g1,g2);
			
	elif len(gene2) == 1:
		gene2=gene2[0];
		lastWord = gene2[-1];
		if lastWord == "*":
			data = service.listExistsGenesRegex(sampleid,gene2[0:-1]);
			geneCount = len(data);
			if geneCount == 1:
				gene = data[0];
				data = service.getExprdataByGene(sampleid,data[0]);
				
		else:
			data = service.getExprdataByGene(sampleid,gene2);
			gene=gene2;
			if data== None:
				geneCount=0;
			else:
				geneCount=1;
	maxval=0;
	if geneCount==1:
		for i in data:
			if data[i] > maxval:
				maxval =data[i];

	return JsonResponse({"res":data,"count":geneCount,"gene":gene,"maxval":maxval});

def getClusterCellids(request):
	sampleid = request.POST.get("sampleid");
	return JsonResponse({"res":sampleid});



def savecluster(request):
	sampleid = request.POST.get("sampleid");
	name = request.POST.get("name");
	ctype = request.POST.get("type");
	comment = request.POST.get("comment");
	cells = request.POST.get("cells");
	cells = cells.split(",");
	cells2 = [];
	for i in cells:
		cells2.append(int(i));

	marks = request.POST.get("marks");
	marks = marks.split(",");
	marks2 = [];
	for i in marks:
		marks2.append( i );


	negmarks = request.POST.get("negmarks");
	negmarks = negmarks.split(",");
	negmarks2 = [];
	for i in negmarks:
		negmarks2.append( i );

	data = service.savecluster(sampleid,name,ctype,cells2,comment,marks2,negmarks2);

	return JsonResponse(data);





def queryClstrCellsAndLabelByCid(request):

	cid = request.POST.get("cid");
	data = service.queryClstrCellsAndLabelByCid(cid);


	return JsonResponse(data);



def getSampleLists(request):
	userid = request.POST.get("userid");

	samples = service.getAllSampleInfo(userid);
	
	jsonres=dict();
	jsonres["samples"]=samples;

	return JsonResponse(jsonres)



def getClusterClassification(request):

	clstrType= request.POST.get("clstrType");

	data = service.getClusterClassification(clstrType);

	return JsonResponse({"clstrTypes":data})


def updatecluster(request):
	target = request.POST.get("target");
	clstrid = request.POST.get("clstrid");

	if target == "POS":
		x = request.POST.get("x");
		y = request.POST.get("y");
		
		res = service.updateClusterPostition(clstrid,x,y);
	elif target =="NAME":
		newname = request.POST.get("name");
		res = service.updateClusterName(clstrid,newname);
	elif target =="prerender":
		val = request.POST.get("val");
		res = service.updateClusterIsPreRender(clstrid,val);

	elif target =="MARKS":
		val = request.POST.get("marks");
		val = val.split(",")
		res = service.updateClusterMarks(clstrid,val);

	elif target =="NEGMARKS":
		val = request.POST.get("marks");
		val = val.split(",")
		res = service.updateClusterNegMarks(clstrid,val);

	elif target =="BOTHMARKS":
		val = request.POST.get("marks");
		negmarks =request.POST.get("negmarks");
		val = val.split(",")
		negmarks = negmarks.split(",")
		res0= service.updateClusterMarks(clstrid,negmarks);
		res = service.updateClusterNegMarks(clstrid,val);

	elif target =="COLOR":
		val = request.POST.get("val");
		val = val.strip();
		res = service.updateClusterColor(clstrid,val);


	return JsonResponse({"res":res});


def deleteCluster(request):
	clstrid = request.POST.get("clstrid");
	res = service.deleteCluster(clstrid);

	return JsonResponse({"res":res});

def contrast(request):
	cells = request.POST.get("cells");
	target = request.POST.get("target");
	sampleid = request.POST.get("sampleid");

	cells = cells.split(",");

	if target=="ALL":
		#data = service.contrastwithrest(sampleid,cells);
		data = service.doContrast(sampleid,cells,"","contrastwithrest");
	else:
		clstrid = target
		data = service.doContrast(sampleid,cells,clstrid,"contrastCellsVsClstr");
		
		#data = service.contrast()


	return JsonResponse({"res":data});

def getContrastResult(request):
	resultid = request.POST.get("resultid");

	res = service.getContrastResult(resultid);

	return JsonResponse(res);


def contrast2(request):
	sampleid = request.POST.get("sampleid");
	clstr= request.POST.get("clstr");
	target = request.POST.get("target");
	
	if target=="ALL":
		cells = service.getClusterCellsById(clstr);
		data = service.doContrast(sampleid,cells,"","contrastwithrest");
	else:
		cells = service.getClusterCellsById(clstr);
		clstrid = target
		data = service.doContrast(sampleid,cells,clstrid,"contrastCellsVsClstr");

	return JsonResponse({"res":data});


def contrastGeneSearch(request):
	sampleid = request.POST.get("sampleid");
	data1 = request.POST.get("data1");
	data2 = request.POST.get("data2");
	gene =request.POST.get("gene");
	dttype = request.POST.get("dttype");

	if dttype =='cid':
		name1 = service.getClusterNameById(data1);
		data1 = service.getClusterCellsById(data1);
		
	else:
		data1=service.strarrayToIntarray(data1);
		name1 ='Selected Cells'


	if data2 =="ALL":
		
		data2 = service.getClusterRestCells(sampleid,data1);
		name2= "Others"

	else:
		name2 = service.getClusterNameById(data2);
		data2 = service.getClusterCellsById(data2);


	plotdata,expr,maxval = service.contrastGeneSearch(gene,data1,data2,sampleid,name1,name2);

	return JsonResponse({"expr":expr,"plot":plotdata,"maxval":maxval});




def getAllClusterStudies(request):

	data = service.getAllClusterStudies();

	return JsonResponse({"data":data});


def getAllTissueByStudies(request):
	study = request.POST.get("study");
	data = service.getAllTissueByStudies(study);

	return JsonResponse({"data":data});


def getAllClusterTypesByStudyAndTissues(request):
	study = request.POST.get("study");
	tissue = request.POST.get("tissue");
	tissue = tissue.split("//,")
	data = service.getAllClusterTypesByStudyAndTissues(study,tissue);

	return JsonResponse({"data":data});



def queryComparePlotData(request):

	study = request.POST.get("study");
	tissues = request.POST.get("tissues");
	tissues = tissues.split("//,")
	gene = request.POST.get("gene");

	cluster = request.POST.get("cluster");
	clusterType = request.POST.get("clusterType");


	data = service.queryComparePlotData(study,tissues,gene,cluster,clusterType);


	return JsonResponse(  data );



def getCellsByStudyAndTissueAndclstrtypeAndClstr(request,study,tissue,clstrType,clstr):

	data = service.getCellsByStudyAndTissueAndclstrtypeAndClstr(study,tissue,clstrType,clstr);


	return JsonResponse(data);