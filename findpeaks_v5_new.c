#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXLEN 1000

//Global variables
int * traitids, numtraits=0, gcounter=1;
fpos_t *traitstarts;
char ** ratioID, ** substrate, ** product;
int ** positions;   //Matrix Chrom x pos
double ** pvals;    //Matrix Chrom x pos
int * maxindices;
char buf[MAXLEN];

void skipheader(FILE* fp){
  char buf[MAXLEN];
  fgets(buf,MAXLEN,fp);    
}

void read_traitids(FILE * fp){
  int curtrait=-1, traitid;
  char tchar1[100], tchar2[100], tchar3[100]; 
  fpos_t position;	
  
  printf("Reading traits...");fflush(stdout);
  skipheader(fp);  
  numtraits=0;
  fgetpos(fp,&position);
	
  while(fgets(buf,MAXLEN,fp)){
     sscanf(buf,"%*[^\t]\t%*[^\t]\t%*d\t%*d\t%*d\t%*d\t%*d\t%*lf\t%*lf\t%*lf\t%[^\t]\t%[^\t]\t%*[^\t]\t%[^\t]\t%*[^\t]\t%d\n", tchar2,tchar3,tchar1,&traitid);
    if(traitid!=curtrait){
      numtraits++;
      traitids=(int*) realloc(traitids,numtraits*sizeof(int));
      traitids[numtraits-1]=traitid;
      traitstarts=(fpos_t*) realloc(traitstarts,numtraits*sizeof(fpos_t));
      traitstarts[numtraits-1]=position;
      ratioID=(char**) realloc(ratioID,numtraits*sizeof(char*));
      ratioID[numtraits-1]=(char*) malloc((strlen(tchar1)+1)*sizeof(char));  
      strcpy(ratioID[numtraits-1],tchar1);  
      substrate=(char**) realloc(substrate,numtraits*sizeof(char*));       
      substrate[numtraits-1]=(char*) malloc((strlen(tchar2)+1)*sizeof(char));  
      strcpy(substrate[numtraits-1],tchar2);  
      product=(char**) realloc(product,numtraits*sizeof(char*));  
      product[numtraits-1]=(char*) malloc((strlen(tchar3)+1)*sizeof(char));  
      strcpy(product[numtraits-1],tchar3);  
      curtrait=traitid; 
      //printf("Traitid=%d, ratioID=%s, substrate=%s, product=%s\n",traitid,tchar1,tchar2,tchar3);  
    }        
    fgetpos(fp,&position);
  }  
  rewind(fp);
  printf("%d traits found\n",numtraits);fflush(stdout);
  //for(int i=0;i<numtraits;i++){
  //  printf("Trait %d starts at %d\n",traitids[i],(int) traitstarts[i]);
  //   fflush(stdout);	  
  //}
}

void filterpeaks(int t, int traitid, int numchrom,int winsize,double cothreshold,int quorum, FILE * fp){
  int curchrom,i;

  //fprintf(fp,"Trait %d\n",traitid);
  //Filter all peaks
  for(curchrom=0;curchrom<numchrom;curchrom++){
    int * peaks=NULL, curpeaks=0;
    for(i=0;i<maxindices[curchrom];i++){
      if(pvals[curchrom][i] > 4.0){
        int j, start, end, cq;  
        start = (int)(ceil(i-(winsize/2)));
        if(start < 0) start=0;            
        end = (int)(floor(i+(winsize/2)));
        if(end >= maxindices[curchrom]) end=maxindices[curchrom]-1;
        cq=0;
        for(j=start;j<=end;j++){
          if(pvals[curchrom][j] > cothreshold)
            cq++;              
        }
        if(cq >=quorum){
          //Add to current list of peaks
          curpeaks++;
          peaks=(int*) realloc(peaks,curpeaks*sizeof(int));
          peaks[curpeaks-1]=i;  
          
          /*
          printf("%d\t%d\t%g:\n",curchrom+1,positions[curchrom][i],pvals[curchrom][i]);
            for(j=start;j<=end;j++){
            if(pvals[curchrom][j] > cothreshold)
                printf("\t%d\t%g\n",positions[curchrom][j],pvals[curchrom][j]);              
            }            
          */  
            
        }            
      }
    }        
    if(peaks){
      int prevpeak;  
      int j, start, end, start2, end2;  
      
      /*        
      for(i=0;i<curpeaks;i++)
        printf("%d ",peaks[i]);
      printf("\n");
      for(i=0;i<curpeaks;i++)
        printf("%d ",positions[curchrom][peaks[i]]);
      printf("\n");
      */  
        
      //Cluster peaks that are less than 10,000bp apart
      prevpeak=peaks[0];
      start=prevpeak;
        
      for(i=1;i<curpeaks;i++){
        //Check if the next peak is less than 10000bp further
        if( positions[curchrom][peaks[i]]-positions[curchrom][prevpeak] < 10000){
          prevpeak=peaks[i];   
        }
        else{
          //Finish previous region
          end = prevpeak; 
          //Print chrom,traitnr,metID,substr,prod,begin,end,folder
          fprintf(fp,"%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\tfolder\n",gcounter++,curchrom+1,traitid,ratioID[t],substrate[t],product[t],positions[curchrom][start]-10000,positions[curchrom][end]+10000);  
          start2 = (int)(ceil(start-(winsize/2)));
          if(start2 < 0) start2=0;            
          end2 = (int)(floor(start+(winsize/2)));
          if(end2 >= maxindices[curchrom]) end2=maxindices[curchrom]-1;            
          for(j=start2;j<=end2;j++){  
            if(pvals[curchrom][j] > cothreshold)
                fprintf(fp,"\t%d\t%g\n",positions[curchrom][j],pvals[curchrom][j]);              
          }
          
          //Start new one            
          prevpeak=peaks[i];  
          start=prevpeak;
        }
      }  
      end = prevpeak; 
      fprintf(fp,"%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\tfolder\n",gcounter++,curchrom+1,traitid,ratioID[t],substrate[t],product[t],positions[curchrom][start]-10000,positions[curchrom][end]+10000);  
      start2 = (int)(ceil(start-(winsize/2)));
      if(start2 < 0) start2=0;            
      end2 = (int)(floor(start+(winsize/2)));
      if(end2 >= maxindices[curchrom]) end2=maxindices[curchrom]-1;            
      for(j=start2;j<=end2;j++){  
        if(pvals[curchrom][j] > cothreshold)
            fprintf(fp,"\t%d\t%g\n",positions[curchrom][j],pvals[curchrom][j]);              
      }
      free(peaks);
    }
  }    
}

void process_file(FILE * fp,int winsize,double cothreshold,int quorum, char * outfile){
  int i,t,curtrait,curchrom, idx, curpos,curnum;
  double curpval;  
  long traitstart;
  FILE * op;
    
  op=fopen(outfile,"w");  
  //Process each traitid   
  skipheader(fp);    
  for(t=0;t<numtraits;t++){    
    int numchrom=0,tmp=0,traitid;

    printf("Starting trait %d\n",traitids[t]);fflush(stdout);
    curtrait=traitids[t];
	  
    fsetpos(fp,&(traitstarts[t]));
	  
    //Check number of chromosomes  
    while(fgets(buf,MAXLEN,fp)){
      sscanf(buf,"%*[^\t]\t%*[^\t]\t%d\t%*d\t%*d\t%*d\t%*d\t%*lf\t%*lf\t%*lf\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%d\n", &tmp,&traitid);
      if(traitid==curtrait){
        if(tmp > numchrom){
          numchrom=tmp;
        }            
      }
      else{
        //fseek(fp,-(strlen(buf)+1),SEEK_CUR);
        break;
      }       
    }    
    printf("%d chromosomes\n",numchrom);    

    //Check maximum index for each chromosome
    maxindices=(int*) calloc(numchrom,sizeof(int));
    for(i=0;i<numchrom;i++) maxindices[i]=0;
    curchrom=1;
    tmp=0;
    curnum=0;
    //Reset file pointer to start of current trait to analyze
    //fseek(fp,traitstart,SEEK_SET);
    fsetpos(fp,&(traitstarts[t]));
    
    while(fgets(buf,MAXLEN,fp)){
      sscanf(buf,"%*[^\t]\t%*[^\t]\t%d\t%*d\t%*d\t%*d\t%*d\t%*lf\t%*lf\t%*lf\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%d\n", &tmp,&traitid);
      if(traitid==curtrait){
        if(tmp > curchrom){
          maxindices[curchrom-1]=curnum;
          curchrom=tmp;  
          curnum=1;
        }       
        else curnum++;    
      }
      else{
        //fseek(fp,-(strlen(buf)+1),SEEK_CUR);  
        break;
      }    
    }
    maxindices[curchrom-1]=curnum;    
    printf("Maxindices done\n");fflush(stdout);
    for(i=0;i<numchrom;i++){
      printf("Max[%d]=%d\n",i,maxindices[i]);
    }

    //Allocate data matrices
    printf("Allocating data...");fflush(stdout);
    positions=(int**) calloc(numchrom,sizeof(int*));
    pvals=(double**) calloc(numchrom,sizeof(double*));
    for(i=0;i<numchrom;i++){
      positions[i]=(int*) calloc(maxindices[i],sizeof(int));
      pvals[i]=(double*) calloc(maxindices[i],sizeof(double));  
    }      
    printf("done\n");fflush(stdout);   

    //Read in all values
    printf("Reading values....");fflush(stdout);
    //fseek(fp,traitstart, SEEK_SET);
    fsetpos(fp,&(traitstarts[t]));
    
    curchrom=1;
    idx=-1;
    while(fgets(buf,MAXLEN,fp)){
      sscanf(buf,"%*[^\t]\t%*[^\t]\t%d\t%d\t%*d\t%*d\t%*d\t%*lf\t%*lf\t%lf\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%d\n", &tmp,&curpos,&curpval,&traitid);
      if(traitid==curtrait){
        if(tmp > curchrom){
          idx=0;
          curchrom=tmp;
        }       
        else idx++;    
        positions[curchrom-1][idx]=curpos;
        pvals[curchrom-1][idx]=curpval;
        if(curpval > 6.0){
          printf("%d\t%d\t%g\n",curchrom,curpos,curpval);
        }
      }
      else{
        //fseek(fp,-(strlen(buf)+1),SEEK_CUR);  
        break;          
      }
    }      
    printf("done\n");fflush(stdout);    

    //Now process data
    printf("Filtering peaks...");fflush(stdout);
    filterpeaks(t,traitids[t], numchrom,winsize,cothreshold,quorum, op);    
    printf("done\n");fflush(stdout);    

    //Free all memory
    for(i=0;i<numchrom;i++){
      free(positions[i]);
      free(pvals[i]);      
    } 
    free(positions);
    free(pvals);
    free(maxindices);
  }    
  fclose(op);
}

int main(int argc, char * argv[]){
  FILE * fp;
  char outfile[250];  
  int winsize=20, quorum=5, i;
  double cothreshold=3.0;
    
  if(argc < 2){
    printf("Syntax: findpeaks inputfile [winsize] [cothreshold] [quorum]\n");
    exit(1);
  }      
    
  fp=fopen(argv[1],"r");
  if(!fp){
    printf("Could not open %s\n",argv[1]);
    exit(1);
  }
  sprintf(outfile,"%s.regions",argv[1]);
  
  if(argc > 2){
      winsize=atoi(argv[2]);
  }
  if(argc > 3){
      cothreshold=atof(argv[3]);
  }
  if(argc > 4){
      quorum=atoi(argv[4]);
  }
  
  //First read all trait ids
  read_traitids(fp);  
  
  //Read file trait by trait
  process_file(fp,winsize,cothreshold,quorum, outfile);
  
  free(traitids);
  for(i=0;i<numtraits;i++){
     free(ratioID[i]);
     free(substrate[i]);
     free(product[i]);
  }
  free(ratioID);
  free(substrate);
  free(product);
}
