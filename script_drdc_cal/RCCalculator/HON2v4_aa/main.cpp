#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
int indexs[21];
int aa(int,int,int);
double radical(double, double, int, int, int,double*);
double conserve(double, double, int, int, int,double*);
double synonymousNew(double, double, int, int, int,double*);
double nonsynonymousNew(double, double, int, int, int,double*);
int change(char);
int dna(char);
double synonymous(double,int,int,int,double*);
double estimate();
double ComputeQab(int, int, int, int, int, int, double, double, double, double, int*, int*, int*, int*);
int SynonymousChange(int, int, int, int, int, int);
int HasTransVersion(int, int, int, int, int, int);
int IsTransVersion(int, int);
int IsRadical(int, int, int, int, int, int, int*);
int main(int argv, char *argc[])
{
	FILE *fptr1,*fptr2,*fptr3,*fptr4,*fptr5,*fptr10,*fptr20,*fptr21;
	int division,flag,i,j,k,k0,p,q,x1,x2,x3,y1,y2,y3,a1,a2,a3,a4;
	double tempr1, tempc1;
	char ***sequence;
	int ***sequenceN, **sequenceA;
	char **seqName,aChar;
	double *S,*N,**s,**n;
	double *C,*R,**c,**r;
	double ratio,temp,temps,tempn;
	double tempc,tempr;
	double sumw;
	char aminoAcid[20]={
		'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'
	}
	;
	int seqLength,seqNum;
	double **pn,**ps,**SS, **NN;
	double **pr,**pc,**CC, **RR;
	double ***ss, ***nn, **var_ps, **var_pn;
	double ***cc, ***rr, **var_pc, **var_pr;
	double ****cov_n, ****cov_s;
	double ****cov_r, ****cov_c;
	double KaKsRatio, DrDcRatio, Dr, Dc;
	double ** frequency;
	double * frequencyAll;
	double ** frequencyAllCodon;
	int ii, jj, kk, ll, ii1, jj1, kk1; /*loop variables*/
	int diff1;
	int diff2;
	int correction = 0;
	int Synonymousp; int Transversionp; int Radicalp;
	double wt[6];
	if((fptr1=fopen(argc[1],"r"))==NULL) {
		printf("Can't open the file!\n"); exit(1);
	}
	fscanf(fptr1,"%d",&seqNum); fscanf(fptr1,"%d",&seqLength);
	if(seqLength!=seqLength/3*3) {
		printf("Non-coding sequences!\n"); exit(1);
	}
	else seqLength=seqLength/3;
	printf("Use corrected N and S? 0 not use, 1 use\n");
	scanf("%d", &correction);
	printf("Input transition to transvertion ratio (0.5:no Ts/Tv bias )");
	scanf("%lf", &ratio);
	ratio=ratio*2.0;
	printf("Input Ka/Ks rate ratio");
	scanf("%lf", &KaKsRatio);
	/* printf("Use proportional difference (1) or JC distance (2)?");
	scanf("%d", &flag);
	*/
	flag=1;
	printf("Choose amino acid group files:\n");
	printf(" 1. Charge\n");
	printf(" 2. Polarity\n");
	printf(" 3. Miyata and Yasunaga\n");
	printf(" 4. Self-defined\n");
	scanf("%d", &division);
	if (division==1) fptr5=fopen("charge.div","r");
	if (division==2) fptr5=fopen("polarity.div","r");
	if (division==3) fptr5=fopen("MY.div","r");
	if (division==4) fptr5=fopen("self.div","r");
	fscanf(fptr5,"%d",&j);
	indexs[20]=-1;
	for(i=0;i<20;i++) indexs[i]=j;
	for(k=1;k<j;k++) {
		fscanf(fptr5,"%d",&k0);
		fscanf(fptr5,"%c",&aChar);
		for(i=0;i<k0;i++) {
			fscanf(fptr5,"%c",&aChar);
			indexs[change(aChar)]=k;
		}
	}
	fclose(fptr5);
	sequence=(char***)malloc(sizeof(char **)*seqNum);
	for(i=0;i<seqNum;i++) sequence[i]=(char**)malloc(sizeof(char *)*seqLength);
	for(i=0;i<seqNum;i++) {
		for(j=0;j<seqLength;j++) {
			sequence[i][j]=(char*)malloc(sizeof(char)*3);
		}
	}
	sequenceN=(int***)malloc(sizeof(int **)*seqNum);
	for(i=0;i<seqNum;i++) sequenceN[i]=(int**)malloc(sizeof(int *)*seqLength);
	for(i=0;i<seqNum;i++) {
		for(j=0;j<seqLength;j++) {
			sequenceN[i][j]=(int*)malloc(sizeof(int)*3);
		}
	}
	sequenceA=(int**)malloc(sizeof(int *)*seqNum);
	for(i=0;i<seqNum;i++) sequenceA[i]=(int*)malloc(sizeof(int)*seqLength);
	frequency=(double**)malloc(sizeof(double *)*seqNum);
	frequencyAll=(double*)malloc(sizeof(double)*70);
	frequencyAllCodon=(double**)malloc(sizeof(double *)*seqNum);
	for(i=0;i<seqNum;i++) frequency[i]=(double*)malloc(sizeof(double)*70);
	for(i=0;i<seqNum;i++) frequencyAllCodon[i]=(double*)malloc(sizeof(double)*70);
	seqName=(char**)malloc(sizeof(char *)*seqNum);
	for (i=0;i<seqNum;i++) seqName[i]=(char*)malloc(sizeof(char)*20);
	S=(double*)malloc(sizeof(double)*seqNum);
	N=(double*)malloc(sizeof(double)*seqNum);
	C=(double*)malloc(sizeof(double)*seqNum);
	R=(double*)malloc(sizeof(double)*seqNum);
	SS=(double**)malloc(sizeof(double *)*seqNum);
	for(i=0;i<seqNum;i++) SS[i]=(double*)malloc(sizeof(double)*seqLength);
	NN=(double**)malloc(sizeof(double *)*seqNum);
	for(i=0;i<seqNum;i++) NN[i]=(double*)malloc(sizeof(double)*seqLength);
	CC=(double**)malloc(sizeof(double *)*seqNum);
	for(i=0;i<seqNum;i++) CC[i]=(double*)malloc(sizeof(double)*seqLength);
	RR=(double**)malloc(sizeof(double *)*seqNum);
	for(i=0;i<seqNum;i++) RR[i]=(double*)malloc(sizeof(double)*seqLength);
	ps=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) ps[i]=(double*)malloc(sizeof(double)*seqNum);
	pn=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) pn[i]=(double*)malloc(sizeof(double)*seqNum);
	pc=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) pc[i]=(double*)malloc(sizeof(double)*seqNum);
	pr=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) pr[i]=(double*)malloc(sizeof(double)*seqNum);
	s=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) s[i]=(double*)malloc(sizeof(double)*seqNum);
	n=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) n[i]=(double*)malloc(sizeof(double)*seqNum);
	c=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) c[i]=(double*)malloc(sizeof(double)*seqNum);
	r=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) r[i]=(double*)malloc(sizeof(double)*seqNum);
	ss=(double***)malloc(sizeof(double **)*seqNum);
	for (i=0;i<seqNum;i++) ss[i]=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) ss[i][j]=(double*)malloc(sizeof(double)*seqLength);
	}
	nn=(double***)malloc(sizeof(double **)*seqNum);
	for (i=0;i<seqNum;i++) nn[i]=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) nn[i][j]=(double*)malloc(sizeof(double)*seqLength);
	}
	cc=(double***)malloc(sizeof(double **)*seqNum);
	for (i=0;i<seqNum;i++) cc[i]=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) cc[i][j]=(double*)malloc(sizeof(double)*seqLength);
	}
	rr=(double***)malloc(sizeof(double **)*seqNum);
	for (i=0;i<seqNum;i++) rr[i]=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) rr[i][j]=(double*)malloc(sizeof(double)*seqLength);
	}
	cov_n=(double****)malloc(sizeof(double ***)*seqNum);
	for (i=0;i<seqNum;i++) cov_n[i]=(double***)malloc(sizeof(double **)*seqNum);
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) {
			cov_n[i][j]=(double**)malloc(sizeof(double *)*seqNum);
		}
	}
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) {
			for (p=0;p<seqNum;p++) {
				cov_n[i][j][p]=(double*)malloc(sizeof(double)*seqNum);
			}
		}
	}
	cov_s=(double****)malloc(sizeof(double ***)*seqNum);
	for (i=0;i<seqNum;i++) cov_s[i]=(double***)malloc(sizeof(double **)*seqNum);
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) {
			cov_s[i][j]=(double**)malloc(sizeof(double *)*seqNum);
		}
	}
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) {
			for (p=0;p<seqNum;p++) {
				cov_s[i][j][p]=(double*)malloc(sizeof(double)*seqNum);
			}
		}
	}
	cov_r=(double****)malloc(sizeof(double ***)*seqNum);
	for (i=0;i<seqNum;i++) cov_r[i]=(double***)malloc(sizeof(double **)*seqNum);
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) {
			cov_r[i][j]=(double**)malloc(sizeof(double *)*seqNum);
		}
	}
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) {
			for (p=0;p<seqNum;p++) {
				cov_r[i][j][p]=(double*)malloc(sizeof(double)*seqNum);
			}
		}
	}
	cov_c=(double****)malloc(sizeof(double ***)*seqNum);
	for (i=0;i<seqNum;i++) cov_c[i]=(double***)malloc(sizeof(double **)*seqNum);
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) {
			cov_c[i][j]=(double**)malloc(sizeof(double *)*seqNum);
		}
	}
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) {
			for (p=0;p<seqNum;p++) {
				cov_c[i][j][p]=(double*)malloc(sizeof(double)*seqNum);
			}
		}
	}
	var_ps=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) var_ps[i]=(double*)malloc(sizeof(double )*seqNum);
	var_pn=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) var_pn[i]=(double*)malloc(sizeof(double )*seqNum);
	var_pc=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) var_pc[i]=(double*)malloc(sizeof(double )*seqNum);
	var_pr=(double**)malloc(sizeof(double *)*seqNum);
	for (i=0;i<seqNum;i++) var_pr[i]=(double*)malloc(sizeof(double )*seqNum);
	for(i=0;i<seqNum;i++) {
		for(j=0;j<seqNum;j++) {
			pn[i][j]=0.0;
			ps[i][j]=0.0;
			n[i][j]=0.0;
			r[i][j]=0.0;
			pc[i][j]=0.0;
			pr[i][j]=0.0;
			s[i][j]=0.0;
			c[i][j]=0.0;
			r[i][j]=0.0;
			for(p=0;p<seqNum;p++) {
				for(q=0;q<seqNum;q++) {
					cov_s[i][j][p][q]=.0; cov_n[i][j][p][q]=.0;
					cov_c[i][j][p][q]=.0; cov_r[i][j][p][q]=.0;
				}
			}
		}
	}
	for(i=0;i<seqNum;i++) {
		for(j=0;j<seqNum;j++) {
			var_ps[i][j]=.0; var_pn[i][j]=.0;
			s[i][j]=.0; n[i][j]=.0;
			var_pc[i][j]=.0; var_pr[i][j]=.0;
			c[i][j]=.0; r[i][j]=.0;
		}
	}
	for(i=0;i<seqNum;i++) {
		for(j=0;j<seqNum;j++) {
			for(k=0;k<seqLength;k++) {
				ss[i][j][k]=.0; nn[i][j][k]=.0; SS[i][k]=.0; NN[i][k]=.0;
				cc[i][j][k]=.0; rr[i][j][k]=.0; CC[i][j]=.0; RR[i][j]=.0;
			}
		}
	}
	for(i=0;i<seqNum;i++) {
		while(fgetc(fptr1)!='\n'){
		}
		j=0;
		while ((seqName[i][j]=fgetc(fptr1))!='\n'&&j<20) {
			j=j+1;
		}
		if(seqName[i][j]!='\n') {
			seqName[i][j]='\n';while(fgetc(fptr1)!='\n') {
			}
		}
		for(j=0;j<seqLength;j++) {
			sequence[i][j][0]=fgetc(fptr1); sequenceN[i][j][0]=dna(sequence[i][j][0]);
			sequence[i][j][1]=fgetc(fptr1); sequenceN[i][j][1]=dna(sequence[i][j][1]);
			sequence[i][j][2]=fgetc(fptr1); sequenceN[i][j][2]=dna(sequence[i][j][2]);
		}
	}
	fclose(fptr1);
	/*compute the codon frequency*/
	for(i=0;i<seqNum;i++) {
		for(j=0;j<70;j++) {
			frequency[i][j]=.0;
		frequencyAllCodon[i][j] = 0.0;
		}
	}
	for(i=0;i<70;i++) {
		frequencyAll[i]=0.0;
	}
	for(i=0;i<seqNum;i++) {
		for(j=0;j<seqLength;j++) {
			sequenceA[i][j]=aa(sequenceN[i][j][0],sequenceN[i][j][1],sequenceN[i][j][2]);
			frequency[i][aa(sequenceN[i][j][0],sequenceN[i][j][1],sequenceN[i][j][2])]++;
			frequencyAll[aa(sequenceN[i][j][0],sequenceN[i][j][1],sequenceN[i][j][2])]++;
			frequencyAllCodon[i][sequenceN[i][j][0]*16+4*sequenceN[i][j][1]+sequenceN[i][j][2]]++;
		}
	}
	for(i=0;i<seqNum;i++) {
		for(j=0;j<70;j++) {
			frequency[i][j]/=seqLength;
		frequencyAllCodon[i][j]/=seqNum;
		}
	}
	//find frequencies for all symbols in the whole input file
	for(j=0;j<70;j++) {
		frequencyAll[j]/=seqNum;
	}

	if(ratio<0.0001) ratio=estimate();

	int round = 0;
	for(i=0;i<seqNum;i++) {
		S[i]=.0; N[i]=.0; C[i]=.0; R[i]=.0;
		for(j=0;j<seqLength;j++) {
			tempr1 = nonsynonymousNew(ratio, KaKsRatio, sequenceN[i][j][0],sequenceN[i][j][1],sequenceN[i][j][2], frequencyAllCodon[i]);
			tempc1 = synonymousNew(ratio, KaKsRatio, sequenceN[i][j][0],sequenceN[i][j][1],sequenceN[i][j][2], frequencyAllCodon[i]);
			sumw = tempr1+tempc1;
			SS[i][j]=3*tempc1/sumw;
			NN[i][j]=3*tempr1/sumw;

			if (correction == 0) 
			{
				SS[i][j]=synonymous(ratio,sequenceN[i][j][0],sequenceN[i][j][1],sequenceN[i][j][2], frequency[i]);
				NN[i][j]=3.0-SS[i][j];
			}

			S[i]=S[i]+SS[i][j];

			/*tempr1=radical(ratio, KaKsRatio, sequenceN[i][j][0],sequenceN[i][j][1],sequenceN[i][j][2], frequencyAll);
			tempc1=conserve(ratio, KaKsRatio, sequenceN[i][j][0],sequenceN[i][j][1],sequenceN[i][j][2], frequencyAll);
			RR[i][j] = NN[i][j]*tempr1/(tempr1+tempc1);*/
			tempr1 = radical(ratio, KaKsRatio, sequenceN[i][j][0],sequenceN[i][j][1],sequenceN[i][j][2], frequency[i]);
			tempc1 = conserve(ratio, KaKsRatio, sequenceN[i][j][0],sequenceN[i][j][1],sequenceN[i][j][2], frequency[i]);
			sumw = tempr1+tempc1;

			CC[i][j]=NN[i][j]*tempc1/sumw;/*weight them according to NN and their weight*/

			RR[i][j]=NN[i][j]*tempr1/sumw;
			R[i]=R[i]+RR[i][j];
		}
		N[i]=3.0*seqLength-S[i];
		C[i]=N[i]-R[i];
	}
	/* computing s, n, c, and r */
	for(i=0;i<seqNum-1;i++) {
		for(j=i+1;j<seqNum;j++) {
			Dr = Dc = 1.0;
			for(k=0;k<seqLength;k++) {
				if(sequenceN[i][k][0]==sequenceN[j][k][0]) a1=0;
				else {
					a1=1;
				}
				if(sequenceN[i][k][1]==sequenceN[j][k][1]) a2=0;
				else {
					a2=1;
				}
				if(sequenceN[i][k][2]==sequenceN[j][k][2]) a3=0;
				else {
					a3=1;
				}
				a4=a1+a2+a3;
				if(a4==1) {
					if(sequenceA[i][k]!=sequenceA[j][k]) {
						if(indexs[sequenceA[i][k]]==indexs[sequenceA[j][k]]) {
							Dc++;
						}
						else {
							Dr++;
						}
					}
				}
			}
			if (Dc == 0) {
				DrDcRatio = 1.0;printf("WARNING : %d %d has Consertive Change=0, set Dr/Dc=1\n", i, j);
			}
			else DrDcRatio = Dr/Dc;
			for(k=0;k<seqLength;k++) {
				round = 0;
				if(sequenceN[i][k][0]==sequenceN[j][k][0]) a1=0;
				else {
					a1=1;
				}
				if(sequenceN[i][k][1]==sequenceN[j][k][1]) a2=0;
				else {
					a2=1;
				}
				if(sequenceN[i][k][2]==sequenceN[j][k][2]) a3=0;
				else {
					a3=1;
				}
				a4=a1+a2+a3;
				if(a4==0) {
					ss[i][j][k]=.0; nn[i][j][k]=.0; cc[i][j][k]=.0; rr[i][j][k]=.0; goto out10;
				}
				if(a4==1) {
					if(sequenceA[i][k]==sequenceA[j][k]) {
						s[i][j]=s[i][j]+1.0; ss[i][j][k]=1.0; nn[i][j][k]=.0;
					}
					else {
						n[i][j]=n[i][j]+1.0; ss[i][j][k]=.0; nn[i][j][k]=1.0;
						if(indexs[sequenceA[i][k]]==indexs[sequenceA[j][k]]) {
							c[i][j]=c[i][j]+1.0; cc[i][j][k]=1.0; rr[i][j][k]=.0;
						}
						else {
							r[i][j]=r[i][j]+1.0; rr[i][j][k]=1.0; cc[i][j][k]=.0;
						}
					}
					goto out10;
				}
				if(a4==2) {
					wt[0] = 0;
					wt[1] = 0;
					if(a1==0 || a2==0) {
						x1=sequenceN[i][k][0];
						x2=sequenceN[i][k][1];
						x3=sequenceN[j][k][2];
						y1=sequenceN[j][k][0];
						y2=sequenceN[j][k][1];
						y3=sequenceN[i][k][2];
					}
					if(a3==0) {
						x1=sequenceN[i][k][0];
						x2=sequenceN[j][k][1];
						x3=sequenceN[i][k][2];
						y1=sequenceN[j][k][0];
						y2=sequenceN[i][k][1];
						y3=sequenceN[j][k][2];
					}
					if(aa(x1,x2,x3)==20) {
						wt[0] = 0; wt[1] = 1; goto start2;
					}
					wt[0] = ComputeQab(sequenceN[i][k][0], sequenceN[i][k][1], sequenceN[i][k][2], 
						x1, x2, x3, frequencyAll[aa(x1,x2,x3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp);

					if(aa(y1,y2,y3)==20) {
						wt[1] = 0; wt[0] = 1; goto start2;
					}
					wt[1] = ComputeQab(x1,x2,x3,
						sequenceN[j][k][0], sequenceN[j][k][1], sequenceN[j][k][2], 
						frequencyAll[aa(sequenceN[j][k][0],sequenceN[j][k][1],sequenceN[j][k][2])], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp);

					temp = wt[0]+wt[1];
					wt[0] /=temp;
					wt[1] /=temp;//weight of the two paths
start2:
					if(a1==0 || a2==0) {
						x1=sequenceN[i][k][0];
						x2=sequenceN[i][k][1];
						x3=sequenceN[j][k][2];
						y1=sequenceN[j][k][0];
						y2=sequenceN[j][k][1];
						y3=sequenceN[i][k][2];
					}
					if(a3==0) {
						x1=sequenceN[i][k][0];
						x2=sequenceN[j][k][1];
						x3=sequenceN[i][k][2];
						y1=sequenceN[j][k][0];
						y2=sequenceN[i][k][1];
						y3=sequenceN[j][k][2];
					}

					temp=.0; temps=.0; tempn=.0; tempc=.0; tempr=.0;
					if(aa(x1,x2,x3)==20) {
						temp=temp+1.0; 
					}
					if(aa(x1,x2,x3)==sequenceA[i][k]) temps=temps+1.0*wt[0];
					else {
						tempn=tempn+1.0*wt[0];
						if(indexs[aa(x1,x2,x3)]==indexs[sequenceA[i][k]]) tempc=tempc+1.0*wt[0];
						else tempr=tempr+1.0*wt[0];
					}
					if(aa(x1,x2,x3)==sequenceA[j][k]) temps=temps+1.0*wt[0];
					else {
						tempn=tempn+1.0*wt[0];
						if(indexs[aa(x1,x2,x3)]==indexs[sequenceA[j][k]]) tempc=tempc+1.0*wt[0];
						else tempr=tempr+1.0*wt[0];
					}
out15:;
					if(aa(y1,y2,y3)==20) {
						temp=temp+1.0; 
					}
					if(aa(y1,y2,y3)==sequenceA[i][k]) temps=temps+1.0*wt[1];
					else {
						tempn=tempn+1.0*wt[1];
						if(indexs[aa(y1,y2,y3)]==indexs[sequenceA[i][k]]) tempc=tempc+1.0*wt[1];
						else tempr=tempr+1.0*wt[1];
					}
					if(aa(y1,y2,y3)==sequenceA[j][k]) temps=temps+1.0*wt[1];
					else {
						tempn=tempn+1.0*wt[1];
						if(indexs[aa(y1,y2,y3)]==indexs[sequenceA[j][k]]) tempc=tempc+1.0*wt[1];
						else tempr=tempr+1.0*wt[1];
					}
out16:;
					s[i][j]=s[i][j]+temps; n[i][j]=n[i][j]+tempn;
					ss[i][j][k]=temps; nn[i][j][k]=tempn;
					c[i][j]=c[i][j]+tempc; r[i][j]=r[i][j]+tempr;
					cc[i][j][k]=tempc; rr[i][j][k]=tempr;
					goto out10;
				}
				temp=.0; temps=.0; tempn=.0; tempc=.0; tempr=.0;
				if(a4==3) {
					wt[0]=wt[1]=wt[2]=wt[3]=wt[4]=wt[5]=0.0;
					x1=sequenceN[i][k][0];
					x2=sequenceN[i][k][1];
					x3=sequenceN[j][k][2];
					y1=sequenceN[i][k][0];
					y2=sequenceN[j][k][1];
					y3=sequenceN[j][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						wt[0]=0;
					}
					else {
						wt[0] = ComputeQab(sequenceN[i][k][0], sequenceN[i][k][1], sequenceN[i][k][2], 
							x1, x2, x3, frequencyAll[aa(x1,x2,x3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(x1,x2,x3, 
							y1, y2, y3, frequencyAll[aa(y1,y2,y3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(y1,y2,y3, 
							sequenceN[j][k][0], sequenceN[j][k][1], sequenceN[j][k][2], frequencyAll[aa(sequenceN[j][k][0],sequenceN[j][k][1],sequenceN[j][k][2])], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp);
					}
					x1=sequenceN[i][k][0];
					x2=sequenceN[i][k][1];
					x3=sequenceN[j][k][2];
					y1=sequenceN[j][k][0];
					y2=sequenceN[i][k][1];
					y3=sequenceN[j][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						wt[1]=0;
					}
					else {
						wt[1] = ComputeQab(sequenceN[i][k][0], sequenceN[i][k][1], sequenceN[i][k][2], 
							x1, x2, x3, frequencyAll[aa(x1,x2,x3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(x1,x2,x3, 
							y1, y2, y3, frequencyAll[aa(y1,y2,y3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(y1,y2,y3, 
							sequenceN[j][k][0], sequenceN[j][k][1], sequenceN[j][k][2], frequencyAll[aa(sequenceN[j][k][0],sequenceN[j][k][1],sequenceN[j][k][2])], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp);
					}
					x1=sequenceN[i][k][0];
					x2=sequenceN[j][k][1];
					x3=sequenceN[i][k][2];
					y1=sequenceN[i][k][0];
					y2=sequenceN[j][k][1];
					y3=sequenceN[j][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						wt[2]=0;
					}
					else {
						wt[2] = ComputeQab(sequenceN[i][k][0], sequenceN[i][k][1], sequenceN[i][k][2], 
							x1, x2, x3, frequencyAll[aa(x1,x2,x3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(x1,x2,x3, 
							y1, y2, y3, frequencyAll[aa(y1,y2,y3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(y1,y2,y3, 
							sequenceN[j][k][0], sequenceN[j][k][1], sequenceN[j][k][2], frequencyAll[aa(sequenceN[j][k][0],sequenceN[j][k][1],sequenceN[j][k][2])], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp);
					}
					x1=sequenceN[i][k][0];
					x2=sequenceN[j][k][1];
					x3=sequenceN[i][k][2];
					y1=sequenceN[j][k][0];
					y2=sequenceN[j][k][1];
					y3=sequenceN[i][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						wt[3]=0;
					}
					else {
						wt[3] = ComputeQab(sequenceN[i][k][0], sequenceN[i][k][1], sequenceN[i][k][2], 
							x1, x2, x3, frequencyAll[aa(x1,x2,x3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(x1,x2,x3, 
							y1, y2, y3, frequencyAll[aa(y1,y2,y3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(y1,y2,y3, 
							sequenceN[j][k][0], sequenceN[j][k][1], sequenceN[j][k][2], frequencyAll[aa(sequenceN[j][k][0],sequenceN[j][k][1],sequenceN[j][k][2])], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp);

					}
					x1=sequenceN[j][k][0];
					x2=sequenceN[i][k][1];
					x3=sequenceN[i][k][2];
					y1=sequenceN[j][k][0];
					y2=sequenceN[i][k][1];
					y3=sequenceN[j][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						wt[4]=0;
					}
					else{
						wt[4] = ComputeQab(sequenceN[i][k][0], sequenceN[i][k][1], sequenceN[i][k][2], 
							x1, x2, x3, frequencyAll[aa(x1,x2,x3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(x1,x2,x3, 
							y1, y2, y3, frequencyAll[aa(y1,y2,y3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(y1,y2,y3, 
							sequenceN[j][k][0], sequenceN[j][k][1], sequenceN[j][k][2], frequencyAll[aa(sequenceN[j][k][0],sequenceN[j][k][1],sequenceN[j][k][2])], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp);

					}
					x1=sequenceN[j][k][0];
					x2=sequenceN[i][k][1];
					x3=sequenceN[i][k][2];
					y1=sequenceN[j][k][0];
					y2=sequenceN[j][k][1];
					y3=sequenceN[i][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						wt[5]=0;
					}
					else {
						wt[5] = ComputeQab(sequenceN[i][k][0], sequenceN[i][k][1], sequenceN[i][k][2], 
							x1, x2, x3, frequencyAll[aa(x1,x2,x3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(x1,x2,x3, 
							y1, y2, y3, frequencyAll[aa(y1,y2,y3)], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp)*
							ComputeQab(y1,y2,y3, 
							sequenceN[j][k][0], sequenceN[j][k][1], sequenceN[j][k][2], frequencyAll[aa(sequenceN[j][k][0],sequenceN[j][k][1],sequenceN[j][k][2])], ratio, KaKsRatio, DrDcRatio, indexs, &Synonymousp, &Transversionp, &Radicalp);

					}
					tempn=wt[0]+wt[1]+wt[2]+wt[3]+wt[4]+wt[5];
					wt[0]/=tempn;
					wt[1]/=tempn;
					wt[2]/=tempn;
					wt[3]/=tempn;
					wt[4]/=tempn;
					wt[5]/=tempn;
					tempn=0.0;
start6:
					x1=sequenceN[i][k][0];
					x2=sequenceN[i][k][1];
					x3=sequenceN[j][k][2];
					y1=sequenceN[i][k][0];
					y2=sequenceN[j][k][1];
					y3=sequenceN[j][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						temp=temp+1; goto out20; 
					}
					if(aa(x1,x2,x3)==sequenceA[i][k]) temps=temps+1.0*wt[0];
					else {
						tempn=tempn+1.0*wt[0];
						if(indexs[aa(x1,x2,x3)]==indexs[sequenceA[i][k]]) tempc=tempc+1.0*wt[0];
						else tempr=tempr+1.0*wt[0];
					}
					if(aa(y1,y2,y3)==sequenceA[j][k]) temps=temps+1.0*wt[0];
					else {
						tempn=tempn+1.0*wt[0];
						if(indexs[aa(y1,y2,y3)]==indexs[sequenceA[j][k]]) tempc=tempc+1.0*wt[0];
						else tempr=tempr+1.0*wt[0];
					}
					if(aa(x1,x2,x3)==aa(y1,y2,y3)) temps=temps+1.0*wt[0];
					else {
						tempn=tempn+1.0*wt[0];
						if(indexs[aa(x1,x2,x3)]==indexs[aa(y1,y2,y3)]) tempc=tempc+1.0*wt[0];
						else tempr=tempr+1.0*wt[0];
					}
out20:;
					x1=sequenceN[i][k][0];
					x2=sequenceN[i][k][1];
					x3=sequenceN[j][k][2];
					y1=sequenceN[j][k][0];
					y2=sequenceN[i][k][1];
					y3=sequenceN[j][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						temp=temp+1; goto out21;
					}
					if(aa(x1,x2,x3)==sequenceA[i][k]) temps=temps+1.0*wt[1];
					else {
						tempn=tempn+1.0*wt[1];
						if(indexs[aa(x1,x2,x3)]==indexs[sequenceA[i][k]]) tempc=tempc+1.0*wt[1];
						else tempr=tempr+1.0*wt[1];
					}
					if(aa(y1,y2,y3)==sequenceA[j][k]) temps=temps+1.0*wt[1];
					else {
						tempn=tempn+1.0*wt[1];
						if(indexs[aa(y1,y2,y3)]==indexs[sequenceA[j][k]]) tempc=tempc+1.0*wt[1];
						else tempr=tempr+1.0*wt[1];
					}
					if(aa(x1,x2,x3)==aa(y1,y2,y3)) temps=temps+1.0*wt[1];
					else {
						tempn=tempn+1.0*wt[1];
						if(indexs[aa(x1,x2,x3)]==indexs[aa(y1,y2,y3)]) tempc=tempc+1.0*wt[1];
						else tempr=tempr+1.0*wt[1];
					}
out21:;
					x1=sequenceN[i][k][0];
					x2=sequenceN[j][k][1];
					x3=sequenceN[i][k][2];
					y1=sequenceN[i][k][0];
					y2=sequenceN[j][k][1];
					y3=sequenceN[j][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						temp=temp+1; goto out22;
					}
					if(aa(x1,x2,x3)==sequenceA[i][k]) temps=temps+1.0*wt[2];
					else {
						tempn=tempn+1.0*wt[2];
						if(indexs[aa(x1,x2,x3)]==indexs[sequenceA[i][k]]) tempc=tempc+1.0*wt[2];
						else tempr=tempr+1.0*wt[2];
					}
					if(aa(y1,y2,y3)==sequenceA[j][k]) temps=temps+1.0*wt[2];
					else {
						tempn=tempn+1.0*wt[2];
						if(indexs[aa(y1,y2,y3)]==indexs[sequenceA[j][k]]) tempc=tempc+1.0*wt[2];
						else tempr=tempr+1.0*wt[2];
					}
					if(aa(x1,x2,x3)==aa(y1,y2,y3)) temps=temps+1.0*wt[2];
					else {
						tempn=tempn+1.0*wt[2];
						if(indexs[aa(x1,x2,x3)]==indexs[aa(y1,y2,y3)]) tempc=tempc+1.0*wt[2];
						else tempr=tempr+1.0*wt[2];
					}
out22:;
					x1=sequenceN[i][k][0];
					x2=sequenceN[j][k][1];
					x3=sequenceN[i][k][2];
					y1=sequenceN[j][k][0];
					y2=sequenceN[j][k][1];
					y3=sequenceN[i][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						temp=temp+1; 
						goto out23;
					}
					if(aa(x1,x2,x3)==sequenceA[i][k]) temps=temps+1.0*wt[3];
					else {
						tempn=tempn+1.0*wt[3];
						if(indexs[aa(x1,x2,x3)]==indexs[sequenceA[i][k]]) tempc=tempc+1.0*wt[3];
						else tempr=tempr+1.0*wt[3];
					}
					if(aa(y1,y2,y3)==sequenceA[j][k]) temps=temps+1.0*wt[3];
					else {
						tempn=tempn+1.0*wt[3];
						if(indexs[aa(y1,y2,y3)]==indexs[sequenceA[j][k]]) tempc=tempc+1.0*wt[3];
						else tempr=tempr+1.0*wt[3];
					}
					if(aa(x1,x2,x3)==aa(y1,y2,y3)) temps=temps+1.0*wt[3];
					else {
						tempn=tempn+1.0*wt[3];
						if(indexs[aa(x1,x2,x3)]==indexs[aa(y1,y2,y3)]) tempc=tempc+1.0*wt[3];
						else tempr=tempr+1.0*wt[3];
					}
out23:;
					x1=sequenceN[j][k][0];
					x2=sequenceN[i][k][1];
					x3=sequenceN[i][k][2];
					y1=sequenceN[j][k][0];
					y2=sequenceN[i][k][1];
					y3=sequenceN[j][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						temp=temp+1; goto out24;
					}
					if(aa(x1,x2,x3)==sequenceA[i][k]) temps=temps+1.0*wt[4];
					else {
						tempn=tempn+1.0*wt[4];
						if(indexs[aa(x1,x2,x3)]==indexs[sequenceA[i][k]]) tempc=tempc+1.0*wt[4];
						else tempr=tempr+1.0*wt[4];
					}
					if(aa(y1,y2,y3)==sequenceA[j][k]) temps=temps+1.0*wt[4];
					else {
						tempn=tempn+1.0*wt[4];
						if(indexs[aa(y1,y2,y3)]==indexs[sequenceA[j][k]]) tempc=tempc+1.0*wt[4];
						else tempr=tempr+1.0*wt[4];
					}
					if(aa(x1,x2,x3)==aa(y1,y2,y3)) temps=temps+1.0*wt[4];
					else {
						tempn=tempn+1.0*wt[4];
						if(indexs[aa(x1,x2,x3)]==indexs[aa(y1,y2,y3)]) tempc=tempc+1.0*wt[4];
						else tempr=tempr+1.0*wt[4];
					}
out24:;
					x1=sequenceN[j][k][0];
					x2=sequenceN[i][k][1];
					x3=sequenceN[i][k][2];
					y1=sequenceN[j][k][0];
					y2=sequenceN[j][k][1];
					y3=sequenceN[i][k][2];
					if(aa(x1,x2,x3)==20||aa(y1,y2,y3)==20) {
						temp=temp+1; goto out25;
					}
					if(aa(x1,x2,x3)==sequenceA[i][k]) temps=temps+1.0*wt[5];
					else {
						tempn=tempn+1.0*wt[5];
						if(indexs[aa(x1,x2,x3)]==indexs[sequenceA[i][k]]) tempc=tempc+1.0*wt[5];
						else tempr=tempr+1.0*wt[5];
					}
					if(aa(y1,y2,y3)==sequenceA[j][k]) temps=temps+1.0*wt[5];
					else {
						tempn=tempn+1.0*wt[5];
						if(indexs[aa(y1,y2,y3)]==indexs[sequenceA[j][k]]) tempc=tempc+1.0*wt[5];
						else tempr=tempr+1.0*wt[5];
					}
					if(aa(x1,x2,x3)==aa(y1,y2,y3)) temps=temps+1.0*wt[5];
					else {
						tempn=tempn+1.0*wt[5];
						if(indexs[aa(x1,x2,x3)]==indexs[aa(y1,y2,y3)]) tempc=tempc+1.0*wt[5];
						else tempr=tempr+1.0*wt[5];
					}
out25:;
					s[i][j]=s[i][j]+temps; n[i][j]=n[i][j]+tempn;
					ss[i][j][k]=temps; nn[i][j][k]=tempn;
					c[i][j]=c[i][j]+tempc; r[i][j]=r[i][j]+tempr;
					cc[i][j][k]=tempc; rr[i][j][k]=tempr;
					printf("c=%f r=%f, n=%f s=%f\n", tempc, tempr, tempn, temps);
				}
out10:;
			}
		}
	}
	/* computing ps, pn, ds, dn, and cov matrix */
	for(i=0;i<seqNum-1;i++) {
		for(j=i+1;j<seqNum;j++) {
			pn[i][j]=2.0*n[i][j]/(N[i]+N[j]); ps[i][j]=2.0*s[i][j]/(S[i]+S[j]);
			pr[i][j]=2.0*r[i][j]/(R[i]+R[j]); pc[i][j]=2.0*c[i][j]/(C[i]+C[j]);
		}
	}
	for(i=0;i<seqNum-1;i++) {
		for(j=i+1;j<seqNum;j++) {
			pn[j][i]=pn[i][j]; ps[j][i]=ps[i][j];
			pr[j][i]=pr[i][j]; pc[j][i]=pc[i][j];
		}
	}
	for(i=0;i<seqNum-1;i++) {
		for(j=i+1;j<seqNum;j++) {
			for(k=0;k<seqLength;k++) {
				var_ps[i][j]=var_ps[i][j]+4.0*(ss[i][j][k]-.5*ps[i][j]*(SS[i][k]+SS[j][k]))*(ss[i][j][k]-.5*ps[i][j]*(SS[i][k]+SS[j][k]))/(S[i]+S[j])/(S[i]+S[j]);
				var_pn[i][j]=var_pn[i][j]+4.0*(nn[i][j][k]-.5*pn[i][j]*(NN[i][k]+NN[j][k]))*(nn[i][j][k]-.5*pn[i][j]*(NN[i][k]+NN[j][k]))/(N[i]+N[j])/(N[i]+N[j]);
				var_pc[i][j]=var_pc[i][j]+4.0*(cc[i][j][k]-.5*pc[i][j]*(CC[i][k]+CC[j][k]))*(cc[i][j][k]-.5*pc[i][j]*(CC[i][k]+CC[j][k]))/(C[i]+C[j])/(C[i]+C[j]);
				var_pr[i][j]=var_pr[i][j]+4.0*(rr[i][j][k]-.5*pr[i][j]*(RR[i][k]+RR[j][k]))*(rr[i][j][k]-.5*pr[i][j]*(RR[i][k]+RR[j][k]))/(R[i]+R[j])/(R[i]+R[j]);
			}
		}
	}
	for(i=0;i<seqNum-1;i++) {
		for(j=i+1;j<seqNum;j++) {
			var_ps[j][i]=var_ps[i][j]; var_pn[j][i]=var_pn[i][j];
			var_pc[j][i]=var_pc[i][j]; var_pr[j][i]=var_pr[i][j];
		}
	}
	for(i=0;i<seqNum-2;i++) {
		for(j=i+1;j<seqNum;j++) {
			for(p=i;p<seqNum-1;p++) {
				if(p==i) x1=j;
				else x1=p;
				for(q=x1+1;q<seqNum;q++) {
					for(k=0;k<seqLength;k++) {
						cov_s[i][j][p][q]=cov_s[i][j][p][q]+(ss[i][j][k]-ps[i][j]*.5*(SS[i][k]+SS[j][k]))*(ss[p][q][k]-ps[p][q]*.5*(SS[p][k]+SS[q][k]))/((S[i]+S[j])*(S[p]+S[q])/4.0);
						cov_c[i][j][p][q]=cov_c[i][j][p][q]+(cc[i][j][k]-pc[i][j]*.5*(CC[i][k]+CC[j][k]))*(cc[p][q][k]-pc[p][q]*.5*(CC[p][k]+CC[q][k]))/((C[i]+C[j])*(C[p]+C[q])/4.0);
					}
				}
			}
		}
	}
	for(i=0;i<seqNum-2;i++) {
		for(j=i+1;j<seqNum;j++) {
			for(p=i;p<seqNum-1;p++) {
				if(p==i) x1=j;
				else x1=p;
				for(q=x1+1;q<seqNum;q++) {
					for(k=0;k<seqLength;k++) {
						cov_n[i][j][p][q]=cov_n[i][j][p][q]+(nn[i][j][k]-pn[i][j]*.5*(NN[i][k]+NN[j][k]))*(nn[p][q][k]-pn[p][q]*.5*(NN[p][k]+NN[q][k]))/((N[i]+N[j])*(N[p]+N[q])/4.0);
						cov_r[i][j][p][q]=cov_r[i][j][p][q]+(rr[i][j][k]-pr[i][j]*.5*(RR[i][k]+RR[j][k]))*(rr[p][q][k]-pr[p][q]*.5*(RR[p][k]+RR[q][k]))/((R[i]+R[j])*(R[p]+R[q])/4.0);
					}
				}
			}
		}
	}
	for(i=0;i<seqNum-1;i++) {
		for(j=i+1;j<seqNum;j++) {
			cov_s[i][j][i][j]=var_ps[i][j]; cov_n[i][j][i][j]=var_pn[i][j];
			cov_c[i][j][i][j]=var_pc[i][j]; cov_r[i][j][i][j]=var_pr[i][j];
		}
	}
	/* output */
	fptr1=fopen("sn.rst","w");
	fptr10=fopen("cr.rst","w");
	fprintf(fptr1,"\nNumbers of nonsynonymous and synonymous sites\n");
	fprintf(fptr10,"\nNumbers of radical and conservative sites\n");
	for(i=0;i<seqNum;i++) {
		fprintf(fptr1,"otu%3d: %8.3f %8.3f\n",i+1,N[i],S[i]);
		fprintf(fptr10,"otu%3d: %8.3f %8.3f\n",i+1,R[i],C[i]);
	}
	for(i=0;i<seqNum-1;i++) {
		for(j=i+1;j<seqNum;j++) {
			s[j][i]=s[i][j]; n[j][i]=n[i][j];
			c[j][i]=c[i][j]; r[j][i]=r[i][j];
		}
	}
	fprintf(fptr1,"\nNumber of nonsynonymous differences");
	fprintf(fptr10,"\nNumber of radical differences");
	for(i=0;i<seqNum;i++) {
		fprintf(fptr1,"\n");fprintf(fptr10,"\n");
		for(j=0;j<seqNum;j++) {
			fprintf(fptr1,"%8.3f ", n[i][j]);
			fprintf(fptr10,"%8.3f ", r[i][j]);
		}
	}
	fprintf(fptr1,"\n");
	fprintf(fptr1,"\nNumber of synonymous differences");
	fprintf(fptr10,"\n");
	fprintf(fptr10,"\nNumber of conservative differences");
	for(i=0;i<seqNum;i++) {
		fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
		for(j=0;j<seqNum;j++) {
			fprintf(fptr1,"%8.3f ", s[i][j]);
			fprintf(fptr10,"%8.3f ", c[i][j]);
		}
	}
	fprintf(fptr1,"\n");
	fprintf(fptr1,"\nProportion of nonsynonymous differences");
	fprintf(fptr10,"\n");
	fprintf(fptr10,"\nProportion of radical differences");
	for(i=0;i<seqNum;i++) {
		fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
		for(j=0;j<seqNum;j++) {
			fprintf(fptr1,"%8.3f ", pn[i][j]);
			fprintf(fptr10,"%8.3f ", pr[i][j]);
		}
	}
	fprintf(fptr1,"\n");
	fprintf(fptr1,"\nProportion of synonymous differences");
	fprintf(fptr10,"\n");
	fprintf(fptr10,"\nProportion of conservative differences");
	for(i=0;i<seqNum;i++) {
		fprintf(fptr1,"\n");fprintf(fptr10,"\n");
		for(j=0;j<seqNum;j++) {
			fprintf(fptr1,"%8.3f ", ps[i][j]);
			fprintf(fptr10,"%8.3f ", pc[i][j]);
		}
	}
	fprintf(fptr1,"\n");
	fprintf(fptr1,"\nCovarince of proportion of nonsynonymous differences");
	fprintf(fptr10,"\n");
	fprintf(fptr10,"\nCovarince of proportion of radical differences");
	for(i=0;i<seqNum-1;i++) {
		fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
		for(j=i+1;j<seqNum;j++) {
			fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
			for(p=i;p<seqNum-1;p++) {
				fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
				if(p==i) k=j;
				else k=p+1;
				for(q=k;q<seqNum;q++) {
					if(cov_n[i][j][p][q]<.0) cov_n[i][j][p][q]=.0;
					if(cov_r[i][j][p][q]<.0) cov_r[i][j][p][q]=.0;
					fprintf(fptr1,"%d-%d-%d-%d %8.7f ", i,j,p,q,cov_n[i][j][p][q]);
					fprintf(fptr10,"%d-%d-%d-%d %8.7f ", i,j,p,q,cov_r[i][j][p][q]);
				}
			}
		}
	}
	fprintf(fptr1,"\n");
	fprintf(fptr1,"\nCovariance of proportion of synonymous differences");
	fprintf(fptr10,"\n");
	fprintf(fptr10,"\nCovariance of proportion of conservative differences");
	for(i=0;i<seqNum-1;i++) {
		fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
		for(j=i+1;j<seqNum;j++) {
			fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
			for(p=i;p<seqNum-1;p++) {
				fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
				if(p==i) k=j;
				else k=p+1;
				for(q=k;q<seqNum;q++) {
					if(cov_s[i][j][p][q]<.0) cov_s[i][j][p][q]=.0;
					if(cov_c[i][j][p][q]<.0) cov_c[i][j][p][q]=.0;
					fprintf(fptr1,"%d-%d-%d-%d %8.7f ", i,j,p,q,cov_s[i][j][p][q]);
					fprintf(fptr10,"%d-%d-%d-%d %8.7f ", i,j,p,q,cov_c[i][j][p][q]);
				}
			}
		}
	}
	fprintf(fptr1,"\n");
	fprintf(fptr1,"\nJC distance of nonsynonymous substitutions");
	fprintf(fptr10,"\n");
	fprintf(fptr10,"\nJC distance of radical substitutions");
	for(i=0;i<seqNum;i++) {
		fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
		for(j=0;j<seqNum;j++) {
			fprintf(fptr1,"%8.3f ", .0-.75*log(1.0-4.0/3.0*pn[i][j]));
			fprintf(fptr10,"%8.3f ", .0-.75*log(1.0-4.0/3.0*pr[i][j]));
		}
	}
	fprintf(fptr1,"\n");
	fprintf(fptr1,"\nJC distance of synonymous substitutions");
	fprintf(fptr10,"\n");
	fprintf(fptr10,"\nJC distance of conservative substitutions");
	for(i=0;i<seqNum;i++) {
		fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
		for(j=0;j<seqNum;j++) {
			fprintf(fptr1,"%8.3f ", .0-.75*log(1.0-4.0/3.0*ps[i][j]));
			fprintf(fptr10,"%8.3f ", .0-.75*log(1.0-4.0/3.0*pc[i][j]));
		}
	}
	fprintf(fptr1,"\n");
	fprintf(fptr1,"\nCovarince of JC distance of nonsynonymous differences");
	fprintf(fptr10,"\n");
	fprintf(fptr10,"\nCovarince of JC distance of radical differences");
	for(i=0;i<seqNum-1;i++) {
		fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
		for(j=i+1;j<seqNum;j++) {
			fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
			for(p=i;p<seqNum-1;p++) {
				fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
				if(p==i) k=j;
				else k=p+1;
				for(q=k;q<seqNum;q++) {
					if(cov_n[i][j][p][q]<.0) cov_n[i][j][p][q]=.0;
					if(cov_r[i][j][p][q]<.0) cov_r[i][j][p][q]=.0;
					fprintf(fptr1,"%d-%d-%d-%d %8.7f ", i,j,p,q,cov_n[i][j][p][q]/(1.0-4.0/3.0*pn[i][j])/(1.0-4.0/3.0*pn[p][q]));
					fprintf(fptr10,"%d-%d-%d-%d %8.7f ", i,j,p,q,cov_r[i][j][p][q]/(1.0-4.0/3.0*pr[i][j])/(1.0-4.0/3.0*pr[p][q]));
				}
			}
		}
	}
	fprintf(fptr1,"\n");
	fprintf(fptr1,"\nCovariance of JC distance of synonymous differences");
	fprintf(fptr10,"\n");
	fprintf(fptr10,"\nCovariance of JC distance of conservative differences");
	for(i=0;i<seqNum-1;i++) {
		fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
		for(j=i+1;j<seqNum;j++) {
			fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
			for(p=i;p<seqNum-1;p++) {
				fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
				if(p==i) k=j;
				else k=p+1;
				for(q=k;q<seqNum;q++) {
					if(cov_s[i][j][p][q]<.0) cov_s[i][j][p][q]=.0;
					if(cov_c[i][j][p][q]<.0) cov_c[i][j][p][q]=.0;
					fprintf(fptr1,"%d-%d-%d-%d %8.7f ", i,j,p,q,cov_s[i][j][p][q]/(1.0-4.0/3.0*ps[i][j])/(1.0-4.0/3.0*ps[p][q]));
					fprintf(fptr10,"%d-%d-%d-%d %8.7f ", i,j,p,q,cov_c[i][j][p][q]/(1.0-4.0/3.0*pc[i][j])/(1.0-4.0/3.0*pc[p][q]));
				}
			}
		}
	}
	fclose(fptr1);
	fclose(fptr10);
	fptr1=fopen("n.dis","w");
	fptr2=fopen("s.dis","w");
	fptr10=fopen("r.dis","w");
	fptr20=fopen("c.dis","w");
	fprintf(fptr1,"%d\n",seqNum);
	fprintf(fptr2,"%d\n",seqNum);
	fprintf(fptr10,"%d\n",seqNum);
	fprintf(fptr20,"%d\n",seqNum);
	if(flag==1) {
		for(i=0;i<seqNum;i++) {
			fprintf(fptr1,"\n"); fprintf(fptr10,"\n");
			for(j=0;j<seqNum;j++) {
				fprintf(fptr1,"%f ", pn[i][j]); fprintf(fptr10,"%f ", pr[i][j]);
			}
		}
		for(i=0;i<seqNum;i++) {
			fprintf(fptr2,"\n"); fprintf(fptr20,"\n");
			for(j=0;j<seqNum;j++) {
				fprintf(fptr2,"%f ", ps[i][j]); fprintf(fptr20,"%f ", pc[i][j]);
			}
		}
		fprintf(fptr1,"\n");
		fprintf(fptr10,"\n");
		for(i=0;i<seqNum-1;i++) {
			for(j=i+1;j<seqNum;j++) {
				for(p=i;p<seqNum-1;p++) {
					if(p==i) k=j;
					else k=p+1;
					for(q=k;q<seqNum;q++) {
						if(cov_n[i][j][p][q]<.0) cov_n[i][j][p][q]=.0;
						if(cov_r[i][j][p][q]<.0) cov_r[i][j][p][q]=.0;
						fprintf(fptr1,"%f ", cov_n[i][j][p][q]);
						fprintf(fptr10,"%f ", cov_r[i][j][p][q]);
					}
				}
			}
		}
		fprintf(fptr2,"\n");
		fprintf(fptr20,"\n");
		for(i=0;i<seqNum-1;i++) {
			for(j=i+1;j<seqNum;j++) {
				for(p=i;p<seqNum-1;p++) {
					if(p==i) k=j;
					else k=p+1;
					for(q=k;q<seqNum;q++) {
						if(cov_s[i][j][p][q]<.0) cov_s[i][j][p][q]=.0;
						if(cov_c[i][j][p][q]<.0) cov_c[i][j][p][q]=.0;
						fprintf(fptr2,"%f ",cov_s[i][j][p][q]);
						fprintf(fptr20,"%f ",cov_c[i][j][p][q]);
					}
				}
			}
		}
	}
	else {
		for(i=0;i<seqNum;i++) {
			fprintf(fptr1,"\n");fprintf(fptr10,"\n");
			for(j=0;j<seqNum;j++) {
				fprintf(fptr1,"%f ", .0-.75*log(1.0-4.0/3.0*pn[i][j]));
				fprintf(fptr10,"%f ", .0-.75*log(1.0-4.0/3.0*pr[i][j]));
			}
		}
		for(i=0;i<seqNum;i++) {
			fprintf(fptr2,"\n"); fprintf(fptr20,"\n");
			for(j=0;j<seqNum;j++) {
				fprintf(fptr2,"%8.3f ", .0-.75*log(1.0-4.0/3.0*ps[i][j]));
				fprintf(fptr20,"%8.3f ", .0-.75*log(1.0-4.0/3.0*pc[i][j]));
			}
		}
		fprintf(fptr1,"\n");
		fprintf(fptr10,"\n");
		for(i=0;i<seqNum-1;i++) {
			for(j=i+1;j<seqNum;j++) {
				for(p=i;p<seqNum-1;p++) {
					if(p==i) k=j;
					else k=p+1;
					for(q=k;q<seqNum;q++) {
						if(cov_n[i][j][p][q]<.0) cov_n[i][j][p][q]=.0;
						if(cov_r[i][j][p][q]<.0) cov_r[i][j][p][q]=.0;
						fprintf(fptr1,"%f ",cov_n[i][j][p][q]/(1.0-4.0/3.0*pn[i][j])/(1.0-4.0/3.0*pn[p][q]));
						fprintf(fptr10,"%f ",cov_r[i][j][p][q]/(1.0-4.0/3.0*pr[i][j])/(1.0-4.0/3.0*pr[p][q]));
					}
				}
			}
		}
		fprintf(fptr2,"\n");
		fprintf(fptr20,"\n");
		for(i=0;i<seqNum-1;i++) {
			for(j=i+1;j<seqNum;j++) {
				for(p=i;p<seqNum-1;p++) {
					if(p==i) k=j;
					else k=p+1;
					for(q=k;q<seqNum;q++) {
						if(cov_s[i][j][p][q]<.0) cov_s[i][j][p][q]=.0;
						if(cov_c[i][j][p][q]<.0) cov_c[i][j][p][q]=.0;
						fprintf(fptr2,"%f ",cov_s[i][j][p][q]/(1.0-4.0/3.0*ps[i][j])/(1.0-4.0/3.0*ps[p][q]));
						fprintf(fptr20,"%f ",cov_c[i][j][p][q]/(1.0-4.0/3.0*pc[i][j])/(1.0-4.0/3.0*pc[p][q]));
					}
				}
			}
		}
	}
	fclose(fptr20);
	fptr21=fopen("outfile","w");
	fprintf(fptr21,"Ts/Tv=%f\n",ratio/2.0);
	if (division==1) fprintf(fptr21,"Charge change");
	if (division==2) fprintf(fptr21,"Polarity change");
	if (division==3) fprintf(fptr21,"Miyata-Yasunaga change");
	if (division==4) fprintf(fptr21,"Self-defined change");
	fprintf(fptr21,"\n\nNumbers of radical (R) and conservative (C) sites\n");
	for(i=0;i<seqNum;i++) fprintf(fptr21,"otu%3d: %8.3f %8.3f\n",i+1,R[i],C[i]);
	fprintf(fptr21,"\nNumbers of radical (r, above diagonal) and conservative (c, below diagonal) differences");
	for(i=0;i<seqNum;i++) {
		fprintf(fptr21,"\n");
		for(j=0;j<seqNum;j++) {
			if(i==j) fprintf(fptr21,"%8.3f ", 0.0);
			if(j>i) fprintf(fptr21,"%8.3f ", r[i][j]);
			if(j<i) fprintf(fptr21,"%8.3f ", c[i][j]);
		}
	}
	fprintf(fptr21,"\n");
	fprintf(fptr21,"\nProportions of radical differences (above diagonal) and their standard errors (below diagonal)");
	for(i=0;i<seqNum;i++) {
		fprintf(fptr21,"\n"); fprintf(fptr10,"\n");
		for(j=0;j<seqNum;j++) {
			if(i==j) fprintf(fptr21,"%8.3f ", 0.0);
			if(j>i) fprintf(fptr21,"%8.3f ", pr[i][j]);
			if(j<i) fprintf(fptr21,"%8.3f ", sqrt(cov_r[j][i][j][i]));
		}
	}
	fprintf(fptr21,"\n");
	fprintf(fptr21,"\nProportions of conservative differences (above diagonal) and their standard errors (below diagonal)");
	for(i=0;i<seqNum;i++) {
		fprintf(fptr21,"\n"); fprintf(fptr10,"\n");
		for(j=0;j<seqNum;j++) {
			if(i==j) fprintf(fptr21,"%8.3f ", 0.0);
			if(j>i) fprintf(fptr21,"%8.3f ", pc[i][j]);
			if(j<i) fprintf(fptr21,"%8.3f ", sqrt(cov_c[j][i][j][i]));
		}
	}
	fclose(fptr1);
	fclose(fptr2);
	fclose(fptr10);
	fclose(fptr21);
	for(i=0;i<seqNum;i++) {
		for(j=0;j<seqLength;j++) {
			free(sequence[i][j]);
		}
		free(sequence[i]);
	}
	free(sequence);
	for(i=0;i<seqNum;i++) {
		for(j=0;j<seqLength;j++) {
			free(sequenceN[i][j]);
		}
		free(sequenceN[i]);
	}
	free(sequenceN);
	for(i=0;i<seqNum;i++) {
		free(sequenceA[i]);
		free(frequency[i]);
		free(seqName[i]);
		free(SS[i]);free(NN[i]);free(CC[i]);free(RR[i]);
		free(ps[i]);free(pn[i]);free(pc[i]);free(pr[i]);
		free(s[i]);free(n[i]);free(c[i]);free(r[i]);
		free(var_ps[i]);free(var_pn[i]);free(var_pc[i]);free(var_pr[i]);

	}
	free(sequenceA); free(frequency); free(seqName); free(frequencyAll);
	free(SS);free(NN);free(CC);free(RR);free(ps);free(pn);free(pc);free(pr);
	free(s);free(n);free(c);free(r);
	free(var_ps);free(var_pn);free(var_pc);free(var_pr);
	for(i=0;i<seqNum;i++) {
		for(j=0;j<seqNum;j++) {
			free(ss[i][j]);free(nn[i][j]);free(cc[i][j]);free(rr[i][j]);
		}
		free(ss[i]);free(nn[i]);free(cc[i]);free(rr[i]);
	}
	free(ss);free(nn);free(cc);free(rr);
	for (i=0;i<seqNum;i++) {
		for (j=0;j<seqNum;j++) {
			for (p=0;p<seqNum;p++) {
				free(cov_n[i][j][p]);free(cov_s[i][j][p]);
				free(cov_r[i][j][p]);free(cov_c[i][j][p]);
			}
			free(cov_n[i][j]);free(cov_s[i][j]);
			free(cov_r[i][j]);free(cov_c[i][j]);
		}
		free(cov_n[i]);free(cov_s[i]);free(cov_r[i]);free(cov_c[i]);
	}
	free(cov_n);free(cov_s);free(cov_r);free(cov_c);
	free(S);free(C);free(N);free(R);

}
int dna(char a)
{
	if(a=='t'||a=='T') return(0);
	if(a=='c'||a=='C') return(1);
	if(a=='a'||a=='A') return(2);
	if(a=='g'||a=='G') return(3);
}
double synonymous(double ratio, int a, int b, int c, double* frequencyAll)
{
	double x;
	int k;
	k=a*16+b*4+c;
	x=ratio/(ratio+2.0);
	switch (k)
	{
	case 0:
		return(x);
	case 1:
		return(x);
	case 2:
		return(2.0*x);
	case 3:
		return(2.0*x);
	case 4:
		return(1.0);
	case 5:
		return(1.0);
	case 6:
		return(1.0);
	case 7:
		return(1.0);
	case 8:
		return(1.0);
	case 9:
		return(1.0);
	case 10:
		return(.0);
	case 11:
		return(.0);
	case 12:
		return(ratio/(ratio+1.0));
	case 13:
		return(ratio/(ratio+1.0));
	case 14:
		return(.0);
	case 15:
		return(.0);
	case 16:
		return(1.0);
	case 17:
		return(1.0);
	case 18:
		return(1.0+x);
	case 19:
		return(1.0+x);
	case 20:
		return(1.0);
	case 21:
		return(1.0);
	case 22:
		return(1.0);
	case 23:
		return(1.0);
	case 24:
		return(x);
	case 25:
		return(x);
	case 26:
		return(x);
	case 27:
		return(x);
	case 28:
		return(1.0);
	case 29:
		return(1.0);
	case 30:
		return(1.0+1.0/2.0);
	case 31:
		return(1.0+1.0/(ratio+2.0));
	case 32:
		return((ratio+1.0)/(ratio+2.0));
	case 33:
		return((ratio+1.0)/(ratio+2.0));
	case 34:
		return(2.0/(ratio+2.0));
	case 35:
		return(.0);
	case 36:
		return(1.0);
	case 37:
		return(1.0);
	case 38:
		return(1.0);
	case 39:
		return(1.0);
	case 40:
		return(x);
	case 41:
		return(x);
	case 42:
		return(x);
	case 43:
		return(x);
	case 44:
		return(x);
	case 45:
		return(x);
	case 46:
		return(x+1.0/(ratio+1.0));
	case 47:
		return(x+1.0/(ratio+2.0));
	case 48:
		return(1.0);
	case 49:
		return(1.0);
	case 50:
		return(1.0);
	case 51:
		return(1.0);
	case 52:
		return(1.0);
	case 53:
		return(1.0);
	case 54:
		return(1.0);
	case 55:
		return(1.0);
	case 56:
		return(x);
	case 57:
		return(x);
	case 58:
		return(x);
	case 59:
		return(x);
	case 60:
		return(1.0);
	case 61:
		return(1.0);
	case 62:
		return(1.0);
	case 63:
		return(1.0);
	}
}
int aa(int a, int b, int c)
{
	int k;
	k=a*16+b*4+c;
	switch (k)
	{
	case 0:
		return(13);
	case 1:
		return(13);
	case 2:
		return(10);
	case 3:
		return(10);
	case 4:
		return(15);
	case 5:
		return(15);
	case 6:
		return(15);
	case 7:
		return(15);
	case 8:
		return(18);
	case 9:
		return(18);
	case 10:
		return(20);
	case 11:
		return(20);
	case 12:
		return(4);
	case 13:
		return(4);
	case 14:
		return(20);
	case 15:
		return(17);
	case 16:
		return(10);
	case 17:
		return(10);
	case 18:
		return(10);
	case 19:
		return(10);
	case 20:
		return(14);
	case 21:
		return(14);
	case 22:
		return(14);
	case 23:
		return(14);
	case 24:
		return(8);
	case 25:
		return(8);
	case 26:
		return(5);
	case 27:
		return(5);
	case 28:
		return(1);
	case 29:
		return(1);
	case 30:
		return(1);
	case 31:
		return(1);
	case 32:
		return(9);
	case 33:
		return(9);
	case 34:
		return(9);
	case 35:
		return(12);
	case 36:
		return(16);
	case 37:
		return(16);
	case 38:
		return(16);
	case 39:
		return(16);
	case 40:
		return(2);
	case 41:
		return(2);
	case 42:
		return(11);
	case 43:
		return(11);
	case 44:
		return(15);
	case 45:
		return(15);
	case 46:
		return(1);
	case 47:
		return(1);
	case 48:
		return(19);
	case 49:
		return(19);
	case 50:
		return(19);
	case 51:
		return(19);
	case 52:
		return(0);
	case 53:
		return(0);
	case 54:
		return(0);
	case 55:
		return(0);
	case 56:
		return(3);
	case 57:
		return(3);
	case 58:
		return(6);
	case 59:
		return(6);
	case 60:
		return(7);
	case 61:
		return(7);
	case 62:
		return(7);
	case 63:
		return(7);
	default:
		printf("%d %d %d", a, b, c);
		exit(-1);
	}
}
double estimate()
{
	return(1.0);
}
/* CHANGE from amino acid to number*/
int change(char amnio)
{
	if(amnio=='A') return (0);
	if(amnio=='R') return (1);
	if(amnio=='N') return (2);
	if(amnio=='D') return (3);
	if(amnio=='C') return (4);
	if(amnio=='Q') return (5);
	if(amnio=='E') return (6);
	if(amnio=='G') return (7);
	if(amnio=='H') return (8);
	if(amnio=='I') return (9);
	if(amnio=='L') return (10);
	if(amnio=='K') return (11);
	if(amnio=='M') return (12);
	if(amnio=='F') return (13);
	if(amnio=='P') return (14);
	if(amnio=='S') return (15);
	if(amnio=='T') return (16);
	if(amnio=='W') return (17);
	if(amnio=='Y') return (18);
	if(amnio=='V') return (19);
}
double K;
double radical(double ratio, double KaKsRatio, int u0, int v0, int w0, double* frequencyAll)
{
	double x,y,count, x1, y1, tempw;
	double wt[9];
	int wtcount = 0;
	int flag,flag1,i,j,k,k0,u,v,w;
	x=ratio/(ratio+2.0); //weight of codon A
	y=1.0/(ratio+2.0);
	x1 = ratio/(ratio+1.0);
	y1 = 1.0/(ratio+1.0);
	count=.0;
	u=u0; v=v0; w=w0;
	flag=0; flag1=0;

	float K1, K2;
	K2=1; K1=ratio/2;
	K = ratio/2;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u,v0,w0)]*K1;
		else wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u,v0,w0)]*K2;
		wtcount++;
	}

	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u0,v,w0)]*K1;
		else wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u0,v,w0)]*K2;

		wtcount++;
	}

	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u0,v0,w)]*K1;
		else wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u0,v0,w)]*K2;

		wtcount++;
	}

	for (i = 0; i<9;i++) {
		count += wt[i];
	}
	count/=9;
	for (i=0;i<9;i++) {
		wt[i]/=count;
	}
	u=u0; v=v0; w=w0;
	flag=0; flag1=0;

	count = 0.0;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		if(aa(u,v0,w0)==20) {
			flag=flag+1;
			if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) flag1=1;
		}
	}
	u=u0;
	wtcount = 0;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		tempw=wt[wtcount++];
		if((indexs[aa(u,v0,w0)]!=indexs[aa(u0,v0,w0)])&&(indexs[aa(u,v0,w0)]!=-1)) {
			if(flag==0) {
				if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) count=count+x1*tempw;
				else {
					count=count+y*tempw;
				}
			}
		}
	}
	flag=0; flag1=0;
	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		if(aa(u0,v,w0)==20) {
			flag=flag+1;
			if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) flag1=1;
		}
	}
	v=v0;
	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		tempw=wt[wtcount++];

		if((indexs[aa(u0,v,w0)]!=indexs[aa(u0,v0,w0)])&&(indexs[aa(u0,v,w0)]!=-1)) {
			if(flag==0) {
				if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) count=count+x1*tempw;
				else {
					count=count+y1*tempw;
				}
			}
		}
	}
	flag=0; flag1=0;
	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		if(aa(u0,v0,w)==20) {
			flag=flag+1;
			if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) flag1=1;
		}
	}
	w=w0;
	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		tempw=wt[wtcount++];

		if((indexs[aa(u0,v0,w)]!=indexs[aa(u0,v0,w0)])&&(indexs[aa(u0,v0,w)]!=-1)) {
			if(flag==0) {
				if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) count=count+x1*tempw;
				else {
					count=count+y1*tempw;
				}
			}
		}
	}
	return(count);
}


double conserve(double ratio, double KaKsRatio, int u0, int v0, int w0, double* frequencyAll)
{
	double x,y,count, x1, y1, tempw;
	double wt[9];
	int wtcount = 0;
	int flag,flag1,i,j,k,k0,u,v,w;
	float K1, K2;
	x=ratio/(ratio+2.0); //weight of codon A
	y=1.0/(ratio+2.0);
	x1 = ratio/(ratio+1.0);
	y1 = 1.0/(ratio+1.0);
	count=.0;
	u=u0; v=v0; w=w0;

	K2=1;K1=ratio/2;
	flag=0; flag1=0;
	K = ratio/2;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u,v0,w0)]*K1;
		else wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u,v0,w0)]*K2;
		wtcount++;
	}

	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u0,v,w0)]*K1;
		else wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u0,v,w0)]*K2;

		wtcount++;
	}

	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u0,v0,w)]*K1;
		else wt[wtcount] = frequencyAll[aa(u0,v0,w0)]*frequencyAll[aa(u0,v0,w)]*K2;

		wtcount++;
	}

	for (i = 0; i<9;i++) {
		count += wt[i];
	}
	count/=9;
	for (i=0;i<9;i++) {
		wt[i]/=count;
	}

	count = 0.0;
	u=u0; v=v0; w=w0;
	flag=0; flag1=0;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		if(aa(u,v0,w0)==20) {
			flag=flag+1;
			if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) flag1=1;
		}
	}
	u=u0;
	wtcount = 0;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		tempw=wt[wtcount++];
		if(aa(u,v0,w0)!=aa(u0,v0,w0)&&(indexs[aa(u,v0,w0)]==indexs[aa(u0,v0,w0)])&&(indexs[aa(u,v0,w0)]!=-1)) {
			if(flag==0) {
				if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) count=count+x1*tempw;
				else {
					count=count+y*tempw;
				}
			}
		}
	}
	flag=0; flag1=0;
	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		if(aa(u0,v,w0)==20) {
			flag=flag+1;
			if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) flag1=1;
		}
	}
	v=v0;
	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		tempw=wt[wtcount++];

		if(aa(u0,v,w0)!=aa(u0,v0,w0)&&(indexs[aa(u0,v,w0)]==indexs[aa(u0,v0,w0)])&&(indexs[aa(u0,v,w0)]!=-1)) {
			if(flag==0) {
				if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) count=count+x1*tempw;
				else {
					count=count+y1*tempw;
				}
			}
		}
	}
	flag=0; flag1=0;
	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		if(aa(u0,v0,w)==20) {
			flag=flag+1;
			if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) flag1=1;
		}
	}
	w=w0;
	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		tempw=wt[wtcount++];

		if(aa(u0,v0,w)!=aa(u0,v0,w0)&&(indexs[aa(u0,v0,w)]==indexs[aa(u0,v0,w0)])&&(indexs[aa(u0,v0,w)]!=-1)) {
			if(flag==0) {
				if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) count=count+x1*tempw;
				else {
					count=count+y1*tempw;
				}
			}
		}
	}
	return(count);
}

double ComputeQab(int x1, int y1, int z1, int x2, int y2, int z2, double Pi2, double ratio, double w, double r, int * indexs,
	int* Synonymousp, int* Transversionp, int* Radicalp) {
		int Synonymous, Transversion, Radical;
		double rr = 0.0;
		Synonymous = SynonymousChange(x1, y1, z1, x2, y2, z2);
		Transversion = HasTransVersion(x1, y1, z1, x2, y2, z2);
		Radical = IsRadical(x1, y1, z1, x2, y2, z2, indexs);
		*Synonymousp = Synonymous;
		*Transversionp = Transversion;
		*Radicalp = Radical;
		if (Synonymous == 1 && Transversion == 1) rr = Pi2;
		if (Synonymous == 1 && Transversion == 0) rr = ratio;
		if (Synonymous == 0) {
			if (Radical == 0) {
				if (Transversion == 1) rr = w*Pi2/sqrt(r);
				else rr = ratio*w*Pi2/sqrt(r);
			}
			else {
				if (Transversion == 1) rr = w*Pi2*sqrt(r);
				else rr = ratio*w*Pi2*sqrt(r);
			}
		}
		return rr;
}
int SynonymousChange(int x1, int y1, int z1, int x2, int y2, int z2) {
	int aa1, aa2;
	aa1 = aa(x1, y1, z1);
	aa2 = aa(x2, y2, z2);
	if (aa1 == aa2) return 1;
	return 0;
}
int HasTransVersion(int x1, int y1, int z1, int x2, int y2, int z2) {
	if (IsTransVersion(x1, x2) == 1) return 1;
	if (IsTransVersion(y1, y2) == 1) return 1;
	if (IsTransVersion(z1, z2) == 1) return 1;
	return 0;
}
int IsTransVersion(int x1, int x2) {
	switch (x1) {
	case 0: if (x2 == 2 || x2 == 3) return 1;
		break;
	case 1: if (x2 == 2 || x2 == 3) return 1;
		break;
	case 2: if (x2 == 0 || x2 == 1) return 1;
		break;
	case 3: if (x2 == 0 || x2 == 1) return 1;
		break;
	}
	return 0;
}
int IsRadical(int x1, int y1, int z1, int x2, int y2, int z2, int* indexs)
{
	int aa1, aa2;
	/* assume aa1 and aa2 will always be < 20*/
	aa1 = aa(x1, y1, z1);
	aa2 = aa(x2, y2, z2);
	if (indexs[aa1] == indexs[aa2]) return 0;
	return 1;
}





double nonsynonymousNew(double ratio, double KaKsRatio, int u0, int v0, int w0, double* frequencyAll)
{
	double x,y,count, x1, y1, tempw;
	double wt[9];
	int wtcount = 0;
	int flag,flag1,i,j,k,k0,u,v,w;
	float K1, K2;
	x=ratio/(ratio+2.0); //weight of codon A
	y=1.0/(ratio+2.0);
	x1 = ratio/(ratio+1.0);
	y1 = 1.0/(ratio+1.0);
	count=.0;
	u=u0; v=v0; w=w0;

	K2=1;K1=ratio/2;
	flag=0; flag1=0;
	K = ratio/2;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u*16+4*v0+w0]*K1;
		else wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u*16+4*v0+w0]*K2;
		wtcount++;
	}

	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u0*16+4*v+w0]*K1;
		else wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u0*16+4*v+w0]*K2;

		wtcount++;
	}

	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u0*16+4*v0+w]*K1;
		else wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u0*16+4*v0+w]*K2;

		wtcount++;
	}

	for (i = 0; i<9;i++) {
		count += wt[i];
	}
	count/=9;
	for (i=0;i<9;i++) {
		wt[i]/=count;
	}

	u=u0; v=v0; w=w0;
	flag=0; flag1=0;
	count = 0.0;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		if(aa(u,v0,w0)==20) {
			flag=flag+1;
			if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) flag1=1;
		}
	}
	u=u0;
	wtcount = 0;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		tempw=wt[wtcount++];
		if(aa(u,v0,w0)!=aa(u0,v0,w0)) {
			if(flag==0) {
				if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) count=count+x1*tempw;
				else {
					count=count+y*tempw;
				}
			}
		}
	}
	flag=0; flag1=0;
	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		if(aa(u0,v,w0)==20) {
			flag=flag+1;
			if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) flag1=1;
		}
	}
	v=v0;
	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		tempw=wt[wtcount++];

		if(aa(u0,v,w0)!=aa(u0,v0,w0)) {
			if(flag==0) {
				if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) count=count+x1*tempw;
				else {
					count=count+y1*tempw;
				}
			}
		}
	}
	flag=0; flag1=0;
	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		if(aa(u0,v0,w)==20) {
			flag=flag+1;
			if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) flag1=1;
		}
	}
	w=w0;
	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		tempw=wt[wtcount++];

		if(aa(u0,v0,w)!=aa(u0,v0,w0)) {
			if(flag==0) {
				if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) count=count+x1*tempw;
				else {
					count=count+y1*tempw;
				}
			}
		}
	}
	return(count);
}


double synonymousNew(double ratio, double KaKsRatio, int u0, int v0, int w0, double* frequencyAll)
{
	double x,y,count, x1, y1, tempw;
	double wt[9];
	int wtcount = 0;
	int flag,flag1,i,j,k,k0,u,v,w;
	float K1, K2;
	x=ratio/(ratio+2.0); //weight of codon A
	y=1.0/(ratio+2.0);
	x1 = ratio/(ratio+1.0);
	y1 = 1.0/(ratio+1.0);
	count=.0;
	u=u0; v=v0; w=w0;

	K2=1;K1=ratio/2;
	flag=0; flag1=0;
	K = ratio/2;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u*16+4*v0+w0]*K1;
		else wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u*16+4*v0+w0]*K2;
		wtcount++;
	}

	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u0*16+4*v+w0]*K1;
		else wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u0*16+4*v+w0]*K2;

		wtcount++;
	}

	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u0*16+4*v0+w]*K1;
		else wt[wtcount] = frequencyAll[u0*16+4*v0+w0]*frequencyAll[u0*16+4*v0+w]*K2;

		wtcount++;
	}

	for (i = 0; i<9;i++) {
		count += wt[i];
	}
	count/=9;
	for (i=0;i<9;i++) {
		wt[i]/=count;
	}

	count = 0.0;
	u=u0; v=v0; w=w0;
	flag=0; flag1=0;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		if(aa(u,v0,w0)==20) {
			flag=flag+1;
			if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) flag1=1;
		}
	}
	u=u0;
	wtcount = 0;
	for(i=1;i<=3;i++) {
		u=u+1;
		if(u==4) u=0;
		tempw=wt[wtcount++];
		if(aa(u,v0,w0)==aa(u0,v0,w0)) {
			if(flag==0) {
				if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((u0==0&&u==1)||(u0==1&&u==0)||(u0==2&&u==3)||(u0==3&&u==2)) count=count+x1*tempw;
				else {
					count=count+y*tempw;
				}
			}
		}
	}
	flag=0; flag1=0;
	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		if(aa(u0,v,w0)==20) {
			flag=flag+1;
			if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) flag1=1;
		}
	}
	v=v0;
	for(i=1;i<=3;i++) {
		v=v+1;
		if(v==4) v=0;
		tempw=wt[wtcount++];

		if(aa(u0,v,w0)==aa(u0,v0,w0)) {
			if(flag==0) {
				if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((v0==0&&v==1)||(v0==1&&v==0)||(v0==2&&v==3)||(v0==3&&v==2)) count=count+x1*tempw;
				else {
					count=count+y1*tempw;
				}
			}
		}
	}
	flag=0; flag1=0;
	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		if(aa(u0,v0,w)==20) {
			flag=flag+1;
			if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) flag1=1;
		}
	}
	w=w0;
	for(i=1;i<=3;i++) {
		w=w+1;
		if(w==4) w=0;
		tempw=wt[wtcount++];

		if(aa(u0,v0,w)==aa(u0,v0,w0)) {
			if(flag==0) {
				if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) count=count+x*tempw;
				else {
					count=count+y*tempw;
				}
			}
			if(flag==2) {
				count=count+1*tempw;
			}
			if(flag==1&&flag1==1) {
				count=count+.5*tempw;
			}
			if(flag==1&&flag1==0) {
				if((w0==0&&w==1)||(w0==1&&w==0)||(w0==2&&w==3)||(w0==3&&w==2)) count=count+x1*tempw;
				else {
					count=count+y1*tempw;
				}
			}
		}
	}
	return(count);
}
