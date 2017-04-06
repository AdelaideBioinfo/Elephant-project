#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "stdarg.h"
#include <stack>
#include <list>
#include <math.h>
#include <ctime>
#include <cstring>

using namespace std;

typedef struct Interval
{
    string chr;
    int a;
    int b;
    Interval(string ch,int x,int y){
        chr=ch;
        a=x;
        b=y;
    }
} Interval;

string readWord(string s,int& pos){
    char c = s.at(pos);
    string sr="";
    while (pos<s.length() && (c==' ' || c=='\t' || c=='#' || c=='/')) {
        pos++;
        c=s.at(pos);
    }
    while (pos<s.length() && c>=33 && c<=126 && c!='#' && c!='/'){
        sr+=c;
        pos++;
        if (pos<s.length()) c=s.at(pos);
    }
    return sr;
}

int readint(string s,int& pos){
    char c=s.at(pos);
    while (c<'0' || c>'9') {
        pos++;
        c=s.at(pos);
    }
    int l=0;
    while (pos<s.length() && c>='0' && c<='9') {
        l=l*10+(int)(c-'0');
        pos++;
        if (pos<s.length()) c=s.at(pos);
    }
    return l;
}

bool readSign(string s,int& pos){
    string l = readWord(s,pos);
    if (l.compare("+")==0) {
        return true;
    }
    else if (l.compare("-")==0){
        return false;
    }
    cout<<"Not valid sign"<<endl;
    return false;
}

bool intersect(Interval* iv1,Interval* iv2){
    return ((iv1->chr.compare(iv2->chr)==0) && !(iv1->b<(iv2->a) || iv1->a>(iv2->b)));
}

bool contain(list<string> l,string s){
    for (list<string>::iterator iter=l.begin(); iter!=l.end(); iter++) {
        if (s.compare(*iter)==0) {
            return true;
        }
    }
    return false;
}

string merge(string fn,int length,int gap,int dist){
//merge all neighboring intervals that have gap at most 'gap' such that each merged interval has length at most 'length' and finally only intervals have length at least 'dist' are kept.
    ifstream fi(fn.c_str());
    ostringstream g,sd,l;
    l<<length;
    g<<gap;
    sd<<dist;
    string fout=(fn+"_"+l.str()+"_"+g.str()+"_"+sd.str());
    FILE* w = fopen(fout.c_str(),"wt");
    string line;
    list<int> gaplist;
    int pos=0;
    list<Interval*> base;
    getline(fi,line);
    string ch0 = readWord(line,pos);
    int x0 = readint(line,pos);
    int y0 = readint(line,pos);
    Interval* itv = new Interval(ch0,x0,y0);
    base.push_back(itv);
    gaplist.push_back(0);
    int count=0;
    while (getline(fi,line)) {
        pos=0;
        string ch = readWord(line,pos);
        int x = readint(line,pos);
        int y = readint(line,pos);
        itv = new Interval(ch,x,y);
        base.push_back(itv);
        if (ch.compare(ch0)==0) gaplist.push_back(x-y0);
        else gaplist.push_back(0);
        ch0=ch;
        y0=y;
        count++;
    }
    fi.close();
    int* gaparray = new int[count];
    Interval** basearray = new Interval*[count];
    int c=0;
    list<Interval*>::iterator itb = base.begin();
    for (list<int>::iterator iter = gaplist.begin(); iter!=gaplist.end(); iter++) {
        gaparray[c]=*iter;
        basearray[c]=*itb;
        itb++;
        c++;
    }
    int gapindex[count];
    for (int i=0; i<count; i++) {
        gapindex[i]=i;
    }
    //sort gaparray//
    for (int i=0; i<count; i++) {
        for (int j=0; j<count-i-1; j++) {
            if (gaparray[j+1]>0 && gaparray[j]>gaparray[j+1]) {
                int temp = gaparray[j];
                gaparray[j] = gaparray[j+1];
                gaparray[j+1] = temp;
                int tempind = gapindex[j];
                gapindex[j] = gapindex[j+1];
                gapindex[j+1] = tempind;
            }
        }
    }
    //merge the repeats////
    bool* bl = new bool[count];//bl[i]=false means interval i is merged and will not be taken into account in the final result.
    for (int i=0;i<count;i++) bl[i]=true;
    for (int i=0; i<count; i++) {
        int k0 = gapindex[i];
        if (gaparray[k0]>0) {
            int k=k0;
            while (!bl[k]) k++;
            if (basearray[k]->b-basearray[k0-1]->a <=dist || basearray[k]->a-basearray[k0-1]->b<=gap) {
                basearray[k]->a = basearray[k0-1]->a;//merge the read gapindex[i]-1 and gapindex[i]
                bl[k0-1]=false;
            }
        }
    }
    delete[] gaparray;
    for (int i=0; i<count; i++) {
        if (bl[i] && (basearray[i]->b-basearray[i]->a >= length)) {
            fprintf(w,"%s %d %d\n",(basearray[i]->chr).c_str(),basearray[i]->a,basearray[i]->b);
        }
    }
    delete[] bl;
    delete[] basearray;
    fclose(w);
    return fout;
}

void generateBinSeqAll(int nb,string out,string refFile,string* ele,int limit){
//generate binary sequences with reference 'refFile' and interval repeats from 'ele'
//an interval is present iff it has at least 'limit' base pairs within any reference interval
//all repeats are reported, including invariant intervals
    FILE* temp = fopen((out+"_details").c_str(),"wt");
    FILE* phy = fopen((out+".phy").c_str(),"wt");
    FILE* nexus = fopen((out+".nex").c_str(),"wt");
    string line;
    //read all file///
    list<Interval*> base;
    list<string> chromosome_names;
    string old_chnames="";
    ifstream ref(refFile.c_str());
    while (getline(ref,line)) {
        int pos=0;
        string ch = readWord(line,pos);
        if (ch.compare(old_chnames)!=0) {
            chromosome_names.push_back(ch);
            old_chnames=ch;
        }
        int x = readint(line,pos);
        int y = readint(line,pos);
        Interval* iv = new Interval(ch,x,y);
        base.push_back(iv);
    }
    ref.close();
    cout<<"ref "<<base.size()<<endl;
    
    ///get long name
    list<string> names;
    ifstream fname("longname.txt");
    string* id = new string[nb];
    string* longname = new string[nb];
    if (fname.is_open()){
        while (getline(fname,line)) {
            int pos=0;
            names.push_back(readWord(line,pos));
        }
        for (int i=0;i<nb;i++){
            int pos1=ele[i].find_first_of(".");
            id[i]=ele[i].substr(pos1-1,1);
            for (list<string>::iterator iter=names.begin();iter!=names.end();iter++){
                pos1=(*iter).find_first_of("_");
                string ss=(*iter).substr(pos1+1,1);
                if (ss.compare(id[i])==0) longname[i]=*iter;
                else longname[i]=ele[i];
            }
        }
    }
    else{
        cout<<"FILE longname.txt does not exists"<<endl;
        for (int i=0;i<nb;i++){
            longname[i]=ele[i];
        }
    }
    
    int** eleList = new int*[nb];
    for (int i=0;i<nb;i++) eleList[i] = new int[base.size()];
    for (int i=0; i<nb; i++) {
        int nbAb=0;
        ifstream file(ele[i].c_str());
        list<Interval*> ivl;
        while (getline(file,line)){
            int pos=0;
            string ch = readWord(line,pos);
            int x = readint(line,pos);
            int y = readint(line,pos);
            Interval* iv = new Interval(ch,x,y);
            ivl.push_back(iv);;
        }
        list<Interval*>::iterator iter=base.begin();
        list<Interval*>::iterator iter1=ivl.begin();
        int count=0;
        while (iter!=base.end()){
            if (iter1==ivl.end()){
                eleList[i][count]=0;
                nbAb++;
                iter++;
                count++;
            }
            else{
                string currentChr=(*iter)->chr;
                if ((*iter1)->chr.compare(currentChr)==0) {
                    if (intersect(*iter,*iter1)) {
                        int d=0;
                        if ((*iter1)->b <= (*iter)->b){
                            d=(*iter1)->b-max((*iter1)->a,(*iter)->a);
                            iter1++;
                            while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0 && (*iter1)->b<=(*iter)->b) {
                                d+=(*iter1)->b-(*iter1)->a;
                                iter1++;
                            }
                            if (iter1!=ivl.end() && intersect(*iter,*iter1)){
                                d+=(*iter)->b-(*iter1)->a;
                            }
                            eleList[i][count]=d;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                        }
                        else{
                            d=(*iter)->b-max((*iter1)->a,(*iter)->a);
                            eleList[i][count]=d;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                        }
                        if (d<=limit){
                            nbAb++;
                        }
                        count++;
                    }
                    else{
                        while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0 && (*iter1)->b<(*iter)->a) {
                            iter1++;
                        }
                        if (iter1==ivl.end() || !intersect(*iter,*iter1)) {
                            eleList[i][count]=0;
                            nbAb++;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                            count++;
                        }
                    }
                }
                else {
                    while (iter1!=ivl.end() && (!contain(chromosome_names,(*iter1)->chr) )) {
                        iter1++;
                    }
                    while (iter!=base.end() && (iter1==ivl.end() || (*iter)->chr.compare((*iter1)->chr)!=0)){
                        eleList[i][count]=0;
                        nbAb++;
                        iter++;
                        if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                            while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                iter1++;
                            }
                        }
                        count++;
                    }
                }
            }
        }
        file.close();
        cout<<longname[i]<<" "<<ivl.size()<<" "<<nbAb<<endl;
    }
    
    fprintf(nexus,"#NEXUS\nbegin taxa;\ndimensions ntax=%d;\ntaxlabels\n",nb);
    for (int i=0; i<nb; i++) fprintf(nexus,"%s ",longname[i].c_str());
    fprintf(nexus,"\n;\nend;\n\nbegin characters;\ndimensions nchar=");
    fprintf(nexus,"%d;\nformat symbols = \"01\";\nmatrix\n",base.size());
    fprintf(phy,"%d %d\n",nb,base.size());
    //write into result file//
    string* miss = new string[base.size()];
    for (int i=0;i<base.size();i++) miss[i]="";
    string* present = new string[base.size()];
    fprintf(temp,"# chromosome start end miss_taxa present_taxa\n");
    for (int i=0;i<nb;i++){
        int count=0;
        fprintf(nexus,"%s ",longname[i].c_str());
        fprintf(phy,"%s ",longname[i].c_str());
        list<Interval*>::iterator iter=base.begin();
        for (int cc=0;cc<base.size();cc++){
            if (eleList[i][cc]>limit){
                fprintf(nexus,"1");
                fprintf(phy,"1");
                present[cc]+=id[i]+",";
            }
            else{
                fprintf(nexus,"0");
                fprintf(phy,"0");
                miss[cc]+=id[i]+",";
            }
            if (i==nb-1) {
                fprintf(temp,"%d\t%s\t%d\t%d\t%s\t%s\t",count,((*iter)->chr).c_str(),(*iter)->a,(*iter)->b,miss[cc].c_str(),present[cc].c_str());
                fprintf(temp,"\n");
            }
            count++;
            iter++;
        }
        fprintf(nexus,"\n");
    }
    delete[] eleList;
    delete[] id;
    delete[] present;
    delete[] miss;
    fprintf(nexus,";\nend;");
    fclose(nexus);
    fclose(phy);
    fclose(temp);
}

void generateBinSeqVariant(int nb,string out,string refFile,string* ele,int limit){
//generate binary sequences with reference 'refFile' and interval repeats from 'ele'
//an interval is present iff it has at least 'limit' base pairs within any reference interval
//only variant intervals are reported
    FILE* temp = fopen((out+"_details").c_str(),"wt");
    FILE* nexus = fopen((out+".nex").c_str(),"wt");
    FILE* phylip = fopen((out+".phy").c_str(),"wt");
    string line;
    //read all file///
    list<Interval*> base;
    list<string> chromosome_names;
    string old_chnames="";
    ifstream ref(refFile.c_str());
    while (getline(ref,line)) {
        int pos=0;
        string ch = readWord(line,pos);
        if (ch.compare(old_chnames)!=0) {
            chromosome_names.push_back(ch);
            old_chnames=ch;
        }
        int x = readint(line,pos);
        int y = readint(line,pos);
        Interval* iv = new Interval(ch,x,y);
        base.push_back(iv);
    }
    ref.close();
    cout<<"ref "<<base.size()<<endl;
    
    ///get long name
    list<string> names;
    ifstream fname("longname.txt");
    string* id = new string[nb];
    string* longname = new string[nb];
    if (fname.is_open()){
        while (getline(fname,line)) {
            int pos=0;
            names.push_back(readWord(line,pos));
        }
        for (int i=0;i<nb;i++){
            int pos1=ele[i].find_first_of(".");
            id[i]=ele[i].substr(pos1-1,1);
            for (list<string>::iterator iter=names.begin();iter!=names.end();iter++){
                pos1=(*iter).find_first_of("_");
                string ss=(*iter).substr(pos1+1,1);
                if (ss.compare(id[i])==0) longname[i]=*iter;
            }
        }
    }
    else{
        //cout<<"FILE longname.txt does not exists"<<endl;
        for (int i=0;i<nb;i++){
            longname[i]=ele[i];
        }
    }
    
    ////
    
    int** eleList = new int*[nb];
    for (int i=0;i<nb;i++) eleList[i] = new int[base.size()];
    int* index = new int[base.size()];
    for (int i=0; i<nb; i++) {
        int nbAb=0;
        ifstream file(ele[i].c_str());
        list<Interval*> ivl;
        while (getline(file,line)){
            int pos=0;
            string ch = readWord(line,pos);
            int x = readint(line,pos);
            int y = readint(line,pos);
            Interval* iv = new Interval(ch,x,y);
            ivl.push_back(iv);;
        }
        list<Interval*>::iterator iter=base.begin();
        list<Interval*>::iterator iter1=ivl.begin();
        int count=0;
        while (iter!=base.end()){
            if (iter1==ivl.end()){
                eleList[i][count]=0;
                nbAb++;
                if (i==0){
                    index[count]=0;
                }
                iter++;
                count++;
            }
            else{
                string currentChr=(*iter)->chr;
                if ((*iter1)->chr.compare(currentChr)==0) {
                    if (intersect(*iter,*iter1)) {
                        int d=0;
                        if ((*iter1)->b <= (*iter)->b){
                            d=(*iter1)->b-max((*iter1)->a,(*iter)->a);
                            iter1++;
                            while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0 && (*iter1)->b<=(*iter)->b) {
                                d+=(*iter1)->b-(*iter1)->a;
                                iter1++;
                            }
                            if (iter1!=ivl.end() && intersect(*iter,*iter1)){
                                d+=(*iter)->b-(*iter1)->a;
                            }
                            eleList[i][count]=d;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                        }
                        else{
                            d=(*iter)->b-max((*iter1)->a,(*iter)->a);
                            eleList[i][count]=d;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                        }
                        if (i==0){
                            if (d>limit){
                                index[count]=1;
                            }
                            else {
                                index[count]=0;
                            }
                        }
                        if (d<=limit){
                            nbAb++;
                        }
                        count++;
                    }
                    else{
                        while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0 && (*iter1)->b<(*iter)->a) {
                            iter1++;
                        }
                        if (iter1==ivl.end() || !intersect(*iter,*iter1)) {
                            eleList[i][count]=0;
                            nbAb++;
                            if (i==0){
                                index[count]=0;
                            }
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                            count++;
                        }
                    }
                }
                else {
                    while (iter1!=ivl.end() && (!contain(chromosome_names,(*iter1)->chr) )) {
                        iter1++;
                    }
                    while (iter!=base.end() && (iter1==ivl.end() || (*iter)->chr.compare((*iter1)->chr)!=0)){
                        eleList[i][count]=0;
                        nbAb++;
                        if (i==0){
                            index[count]=0;
                        }
                        iter++;
                        if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                            while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                iter1++;
                            }
                        }
                        count++;
                    }
                }
            }
        }
        file.close();
        cout<<longname[i]<<" "<<ivl.size()<<" "<<nbAb<<endl;
    }
    ///compute index//
    for (int i=1;i<nb;i++){
        for (int cc=0;cc<base.size();cc++){
            if ((index[cc]==1 && eleList[i][cc]<=limit) || (index[cc]==0 && eleList[i][cc]>limit) || index[cc]==-1){
                index[cc]=-1;
            }        }
    }
    int ns=0;
    for (int i=0; i<base.size(); i++) {
        if (index[i]==-1) {
            ns++;
        }
    }
    fprintf(phylip,"%d %d\n",nb,ns);
    fprintf(nexus,"#NEXUS\nbegin taxa;\ndimensions ntax=%d;\ntaxlabels\n",nb);
    for (int i=0; i<nb; i++) fprintf(nexus,"%s ",longname[i].c_str());
    fprintf(nexus,"\n;\nend;\n\nbegin characters;\ndimensions nchar=");
    fprintf(nexus,"%d;\nformat symbols = \"01\";\nmatrix\n",ns);
    //fprintf(phylip,"%d %d\n",nb,ns);
    //write into result file//
    string* miss = new string[base.size()];
    for (int i=0;i<base.size();i++) miss[i]="";
    string* present = new string[base.size()];
    fprintf(temp,"# chromosome start end miss_taxa present_taxa\n");
    for (int i=0;i<nb;i++){
        int count=0;
        fprintf(nexus,"%s ",longname[i].c_str());
        fprintf(phylip,"%s ",longname[i].c_str());
        list<Interval*>::iterator iter=base.begin();
        for (int cc=0;cc<base.size();cc++){
            if (index[cc]==-1) {
                if (eleList[i][cc]>limit){
                    fprintf(nexus,"1");
                    fprintf(phylip,"1");
                    present[cc]+=id[i]+",";
                }
                else{
                    fprintf(nexus,"0");
                    fprintf(phylip,"0");
                    miss[cc]+=id[i]+",";
                }
                if (i==nb-1) {
                    fprintf(temp,"%d\t%s\t%d\t%d\t%s\t%s\t",count,((*iter)->chr).c_str(),(*iter)->a,(*iter)->b,miss[cc].c_str(),present[cc].c_str());
                    fprintf(temp,"\n");
                }
                count++;
            }
            iter++;
        };
        fprintf(nexus,"\n");
        fprintf(phylip,"\n");
    }
    delete[] eleList;
    delete[] index;
    delete[] id;
    delete[] miss;
    delete[] present;
    fprintf(nexus,";\nend;");
    fclose(nexus);
    fclose(phylip);
    fclose(temp);
}

void generateBinSeqVariantRestrict(int nb,string out,string refFile,string* ele,int limit){
    //generate binary sequences with reference 'refFile' and interval repeats from 'ele'
    //an interval is present iff it has at least 'limit' base pairs within any reference interval
    //only variant intervals are reported
    //only intervals that are present in the first species provided in 'ele' are reported
    FILE* nexus = fopen((out+".nex").c_str(),"wt");
    FILE* phylip = fopen((out+".phy").c_str(),"wt");
    string line;
    //read all file///
    list<Interval*> base;
    list<string> chromosome_names;
    string old_chnames="";
    ifstream ref(refFile.c_str());
    while (getline(ref,line)) {
        int pos=0;
        string ch = readWord(line,pos);
        if (ch.compare(old_chnames)!=0) {
            chromosome_names.push_back(ch);
            old_chnames=ch;
        }
        int x = readint(line,pos);
        int y = readint(line,pos);
        Interval* iv = new Interval(ch,x,y);
        base.push_back(iv);
    }
    ref.close();
    cout<<"ref "<<base.size()<<endl;
    
    ///get long name
    list<string> names;
    ifstream fname("longname.txt");
    string* id = new string[nb];
    string* longname = new string[nb];
    if (fname.is_open()){
        while (getline(fname,line)) {
            int pos=0;
            names.push_back(readWord(line,pos));
        }
        for (int i=0;i<nb;i++){
            int pos1=ele[i].find_first_of(".");
            id[i]=ele[i].substr(pos1-1,1);
            for (list<string>::iterator iter=names.begin();iter!=names.end();iter++){
                pos1=(*iter).find_first_of("_");
                string ss=(*iter).substr(pos1+1,1);
                if (ss.compare(id[i])==0) longname[i]=*iter;
            }
        }
    }
    else {
     cout<<"FILE longname.txt does not exists"<<endl;
        for (int i=0;i<nb;i++){
            longname[i]=ele[i];
        }
    }
    
    ////
    
    int** eleList = new int*[nb];
    for (int i=0;i<nb;i++) eleList[i] = new int[base.size()];
    int* index = new int[base.size()];
    for (int i=0;i<base.size();i++) index[i]=1;
    for (int i=0; i<nb; i++) {
        int nbAb=0;
        ifstream file(ele[i].c_str());
        list<Interval*> ivl;
        while (getline(file,line)){
            int pos=0;
            string ch = readWord(line,pos);
            int x = readint(line,pos);
            int y = readint(line,pos);
            Interval* iv = new Interval(ch,x,y);
            ivl.push_back(iv);
        }
        list<Interval*>::iterator iter=base.begin();
        list<Interval*>::iterator iter1=ivl.begin();
        int count=0;
        while (iter!=base.end()){
            if (iter1==ivl.end()){
                eleList[i][count]=0;
                if (i==0){
                    index[count]=0;
                }
                if (index[count]!=0) nbAb++;
                iter++;
                count++;
            }
            else{
                string currentChr=(*iter)->chr;
                if ((*iter1)->chr.compare(currentChr)==0) {
                    if (intersect(*iter,*iter1)) {
                        int d=0;
                        if ((*iter1)->b <= (*iter)->b){
                            d=(*iter1)->b-max((*iter1)->a,(*iter)->a);
                            iter1++;
                            while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0 && (*iter1)->b<=(*iter)->b) {
                                d+=(*iter1)->b-(*iter1)->a;
                                iter1++;
                            }
                            if (iter1!=ivl.end() && intersect(*iter,*iter1)){
                                d+=(*iter)->b-(*iter1)->a;
                            }
                            eleList[i][count]=d;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                        }
                        else{
                            d=(*iter)->b-max((*iter1)->a,(*iter)->a);
                            eleList[i][count]=d;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                        }
                        if (i==0 && d<=limit){
                            index[count]=0;
                        }
                        if (d<=limit){
                            if (index[count]!=0) nbAb++;
                        }
                        count++;
                    }
                    else{
                        while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0 && (*iter1)->b<(*iter)->a) {
                            iter1++;
                        }
                        if (iter1==ivl.end() || !intersect(*iter,*iter1)) {
                            eleList[i][count]=0;
                            if (i==0) index[count]=0;
                            if (index[count]!=0) nbAb++;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                            count++;
                        }
                    }
                }
                else {
                    while (iter1!=ivl.end() && (!contain(chromosome_names,(*iter1)->chr) )) {
                        iter1++;
                    }
                    while (iter!=base.end() && (iter1==ivl.end() || (*iter)->chr.compare((*iter1)->chr)!=0)){
                        eleList[i][count]=0;
                        if (i==0) index[count]=0;
                        if (index[count]!=0) nbAb++;
                        iter++;
                        if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                            while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                iter1++;
                            }
                        }
                        count++;
                    }
                }
            }
        }
        file.close();
        cout<<longname[i]<<" "<<ivl.size()<<" "<<nbAb<<endl;
    }
    ///compute index//
    for (int cc=0;cc<base.size();cc++){
        if (index[cc]==1 && eleList[0][cc]<=limit){
            index[cc]=0;
        }
        else if (index[cc]==1 && eleList[0][cc]>limit){
            index[cc]=1;
        }
        else index[cc]=-1;
    }
    for (int i=1;i<nb;i++){
        for (int cc=0;cc<base.size();cc++){
            if ((index[cc]==1 && eleList[i][cc]<=limit) || (index[cc]==0 && eleList[i][cc]>limit) || index[cc]==-2){
                index[cc]=-2;
            }
        }
    }
    int ns=0;
    for (int i=0; i<base.size(); i++) {
        if (index[i]==-2) {
            ns++;
        }
    }
    fprintf(phylip,"%d %d\n",nb,ns);
    fprintf(nexus,"#NEXUS\nbegin taxa;\ndimensions ntax=%d;\ntaxlabels\n",nb);
    for (int i=0; i<nb; i++) fprintf(nexus,"%s ",longname[i].c_str());
    fprintf(nexus,"\n;\nend;\n\nbegin characters;\ndimensions nchar=");
    fprintf(nexus,"%d;\nformat symbols = \"01\";\nmatrix\n",ns);
    //write into result file//
    for (int i=0;i<nb;i++){
        int count=0;
        fprintf(nexus,"%s ",longname[i].c_str());
        fprintf(phylip,"%s ",longname[i].c_str());
        list<Interval*>::iterator iter=base.begin();
        for (int cc=0;cc<base.size();cc++){
            if (index[cc]==-2) {
                if (eleList[i][cc]>limit){
                    fprintf(nexus,"1");
                    fprintf(phylip,"1");
                }
                else{
                    fprintf(nexus,"0");
                    fprintf(phylip,"0");
                }
                count++;
            }
            iter++;
        }
        fprintf(nexus,"\n");
        fprintf(phylip,"\n");
    }
    delete[] eleList;
    delete[] index;
    delete[] id;
    fprintf(nexus,";\nend;");
    fclose(nexus);
    fclose(phylip);
}

void generateBinSeqVariantExcludeCommon(int nb,int rest,string out,string ref,string* ele,int limit){
    //generate binary sequences with reference 'refFile' and interval repeats from 'ele'
    //an interval is present iff it has at least 'limit' base pairs within any reference interval
    //only variant intervals are reported
    //exclude common intervals of the first 'rest' species provided in 'ele'
    FILE* nexus = fopen((out+".nex").c_str(),"wt");
    FILE* phylip = fopen((out+".phy").c_str(),"wt");
    string line;
    //read all file///
    list<Interval*> base;
    list<string> chromosome_names;
    string old_chnames="";
    ifstream allFile(ref.c_str());
    while (getline(allFile,line)) {
        int pos=0;
        string ch = readWord(line,pos);
        if (ch.compare(old_chnames)!=0) {
            chromosome_names.push_back(ch);
            old_chnames=ch;
        }
        int x = readint(line,pos);
        int y = readint(line,pos);
        Interval* iv = new Interval(ch,x,y);
        base.push_back(iv);
    }
    allFile.close();
    
    ///get long name
    list<string> names;
    ifstream fname("longname.txt");
    string* id = new string[nb];
    string* longname = new string[nb];
    if (fname.is_open()){
        while (getline(fname,line)) {
            int pos=0;
            names.push_back(readWord(line,pos));
        }
        for (int i=0;i<nb;i++){
            int pos1=ele[i].find_first_of(".");
            id[i]=ele[i].substr(pos1-1,1);
            for (list<string>::iterator iter=names.begin();iter!=names.end();iter++){
                pos1=(*iter).find_first_of("_");
                string ss=(*iter).substr(pos1+1,1);
                if (ss.compare(id[i])==0) longname[i]=*iter;
            }
        }
    }
    else{
        cout<<"FILE longname.txt does not exists"<<endl;
        for (int i=0;i<nb;i++){
            longname[i]=ele[i];
        }
    }
            ////
    
    int** eleList = new int*[nb];
    for (int i=0;i<nb;i++) eleList[i] = new int[base.size()];
    int* index = new int[base.size()];
    for (int i=0;i<base.size();i++) index[i]=1;
    for (int i=0; i<nb; i++) {
        ifstream file(ele[i].c_str());
        list<Interval*> ivl;
        while (getline(file,line)){
            int pos=0;
            string ch = readWord(line,pos);
            int x = readint(line,pos);
            int y = readint(line,pos);
            Interval* iv = new Interval(ch,x,y);
            ivl.push_back(iv);
        }
        list<Interval*>::iterator iter=base.begin();
        list<Interval*>::iterator iter1=ivl.begin();
        int count=0;
        while (iter!=base.end()){
            if (iter1==ivl.end()){
                eleList[i][count]=0;
                if (i<rest){
                    index[count]=0;
                }
                iter++;
                count++;
            }
            else{
                string currentChr=(*iter)->chr;
                if ((*iter1)->chr.compare(currentChr)==0) {
                    if (intersect(*iter,*iter1)) {
                        int d=0;
                        if ((*iter1)->b <= (*iter)->b){
                            d=(*iter1)->b-max((*iter1)->a,(*iter)->a);
                            iter1++;
                            while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0 && (*iter1)->b<=(*iter)->b) {
                                d+=(*iter1)->b-(*iter1)->a;
                                iter1++;
                            }
                            if (iter1!=ivl.end() && intersect(*iter,*iter1)){
                                d+=(*iter)->b-(*iter1)->a;
                            }
                            eleList[i][count]=d;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                        }
                        else{
                            d=(*iter)->b-max((*iter1)->a,(*iter)->a);
                            eleList[i][count]=d;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                        }
                        if (i<rest && d<=limit){
                            index[count]=0;
                        }
                        count++;
                    }
                    else{
                        while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0 && (*iter1)->b<(*iter)->a) {
                            iter1++;
                        }
                        if (iter1==ivl.end() || !intersect(*iter,*iter1)) {
                            eleList[i][count]=0;
                            if (i<rest) index[count]=0;
                            iter++;
                            if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                                while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                    iter1++;
                                }
                            }
                            count++;
                        }
                    }
                }
                else {
                    while (iter1!=ivl.end() && (!contain(chromosome_names,(*iter1)->chr) )) {
                        iter1++;
                    }
                    while (iter!=base.end() && (iter1==ivl.end() || (*iter)->chr.compare((*iter1)->chr)!=0)){
                        eleList[i][count]=0;
                        if (i<rest) index[count]=0;
                        iter++;
                        if (iter!=base.end() && (*iter)->chr.compare(currentChr)!=0){
                            while (iter1!=ivl.end() && (*iter1)->chr.compare(currentChr)==0) {
                                iter1++;
                            }
                        }
                        count++;
                    }
                }
            }
        }
        file.close();    }
    ///compute index//
    int commonIntervals=0;
    for (int cc=0;cc<base.size();cc++){
        if (index[cc]==0 && eleList[0][cc]<=limit){
            index[cc]=0;
        }
        else if (index[cc]==0 && eleList[0][cc]>limit){
            index[cc]=1;
        }
        else {
            index[cc]=-1;
            commonIntervals++;
        }
    }
    cout<<"common "<<commonIntervals<<endl;
    
    for (int i=1;i<nb;i++){
        int nbAb=0;
        for (int cc=0;cc<base.size();cc++){
            if ((index[cc]==1 && eleList[i][cc]<=limit) || (index[cc]==0 && eleList[i][cc]>limit) || index[cc]==-2){
                index[cc]=-2;
            }
        }
    }
    int ns=0;
    for (int i=0; i<base.size(); i++) {
        if (index[i]==-2) {
            ns++;
        }
    }
    fprintf(phylip,"%d %d\n",nb,ns);
    fprintf(nexus,"#NEXUS\nbegin taxa;\ndimensions ntax=%d;\ntaxlabels\n",nb);
    for (int i=0; i<nb; i++) fprintf(nexus,"%s ",longname[i].c_str());
    fprintf(nexus,"\n;\nend;\n\nbegin characters;\ndimensions nchar=");
    fprintf(nexus,"%d;\nformat symbols = \"01\";\nmatrix\n",ns);
    //write into result file//
    for (int i=0;i<nb;i++){
        int count=0;
        fprintf(nexus,"%s ",longname[i].c_str());
        fprintf(phylip,"%s ",longname[i].c_str());
        list<Interval*>::iterator iter=base.begin();
        int nbAb=0;
        for (int cc=0;cc<base.size();cc++){
            if (index[cc]==-2) {
                if (eleList[i][cc]>limit){
                    fprintf(nexus,"1");
                    fprintf(phylip,"1");
                }
                else{
                    fprintf(nexus,"0");
                    fprintf(phylip,"0");
                    nbAb++;
                }
                count++;
            }
            iter++;
        }
        cout<<longname[i]<<" "<<nbAb<<endl;
        fprintf(nexus,"\n");
        fprintf(phylip,"\n");
    }
    delete[] eleList;
    delete[] index;
    delete[] id;
    fprintf(nexus,";\nend;");
    fclose(nexus);
    fclose(phylip);
}

int main( int argc, char** argv ){
    if (strcmp(argv[1],"merge")==0) {
        int length = atoi(argv[2]);
        int gap = atoi(argv[3]);
        int dist = atoi(argv[4]);
        string ref=argv[5];
        merge(ref,length,gap,dist);
    }
    else if (strcmp(argv[1],"binSeqAll")==0){
        int limit = atoi(argv[2]);
        string fout=argv[3];
        string ref=argv[4];
        string* elem = new string[argc-5];
        for (int i=0; i<argc-5; i++) {
            elem[i]=string(argv[i+5]);
        }
        generateBinSeqAll(argc-5,fout,ref,elem,limit);
        delete[] elem;
    }
    else if (strcmp(argv[1],"binSeqVariant")==0){
        int limit = atoi(argv[2]);
        string fout=argv[3];
        string ref=argv[4];
        string* elem = new string[argc-5];
        for (int i=0; i<argc-5; i++) {
            elem[i]=string(argv[i+5]);
        }
        generateBinSeqVariant(argc-5,fout,ref,elem,limit);
        delete[] elem;
    }
    else if (strcmp(argv[1],"binSeqVariantRestrict")==0){
        int limit = atoi(argv[2]);
        string fout=argv[3];
        string ref=argv[4];
        string* elem = new string[argc-5];
        for (int i=0; i<argc-5; i++) {
            elem[i]=string(argv[i+5]);
        }
        generateBinSeqVariantRestrict(argc-5,fout,ref,elem,limit);
        delete[] elem;
    }
    else if (strcmp(argv[1],"binSeqVariantExcludeCommon")==0){
        int limit = atoi(argv[2]);
        int rest = atoi(argv[3]);
        string fout=argv[4];
        string ref=argv[5];
        string* elem = new string[argc-6];
        for (int i=0; i<argc-6; i++) {
            elem[i]=string(argv[i+6]);
        }
        generateBinSeqVariantExcludeCommon(argc-6,rest,fout,ref,elem,limit);
        delete[] elem;
    }
    else{
        cout<<"unrecognized command "<<argv[1]<<endl;
    }
    return 1;
}

//./binarySequence binSeqVariant 0 output ref.bed A.bed B.bed C.bed D.bed E.bed F.bed
//compute variant intervals of 6 species A, B, .., F


//./binarySequence binSeqVariantRestrict 0 output ref.bed D.bed A.bed B.bed C.bed E.bed F.bed
//compute variant intervals that are present in D of 6 species A, B, ..., F

//./binarySequence binSeqVariantExcludeCommon 0 2 output ref.bed C.bed D.bed A.bed B.bed E.bed
//compute variant intervals of 6 species A, B, ...F that are not present in the common intervals of C and D
