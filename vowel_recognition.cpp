// vowel_recognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include "math.h"
#include "ctype.h"
#include <limits>
#include "stdlib.h"
#include <cstring>

#define N 320 //number of frames in one frame
#define F 5 // stable frames count
#define P 12 

//global variables

/*array for storing the values of ri, ci, allci, reference ci of frame*/
double R[F][P+1], C[F][P+1];
double restoreCi[F][P+1],allCi[50][F][P+1];

int files = 0;// helps to store values of Ci to allCi

/*to store tokhura distances,array for samples*/
double tokhuraDist[5],x[80000];


/* utility function to apply Durbin Algorithm and find The value of ai's
   and find cepstral cofficients and applying raised sine window and 
   storing the ci values to 3D matrix*/
void durbinAlgo()
{
	//applying durbin's algorithm
	double A[F][P+1];
	double sum=0, PI=3.142857142857;
	int f1=0,i1=0,j1=0,i2=0,j2=0;
	double val1,val2,val3,val4;
	double E[13];
	double K[13];
	double Alpha[13][13];
	while(f1<F)
	{
		E[0] = R[f1][0];
		i1=1;
		while(i1<=P)
		{
			sum=0;
			j1=1;
			while(j1<=i1-1)
			{
				val1= Alpha[i1-1][j1] * R[f1][i1-j1];
				sum = sum + val1;	
				++j1;//update j1 by 1
			}
			
			val2 = (R[f1][i1]-sum);
			K[i1]=val2/E[i1-1];
				
			Alpha[i1][i1] = K[i1];
		    j2=1;
			while(j2<=i1-1)
			{
				val3= K[i1]*Alpha[i1-1][i1-j2];
				Alpha[i1][j2] = Alpha[i1-1][j2] - val3;
				++j2;//update j2 by 1
			}
		    val4 = (1-(K[i1]*K[i1]));
			E[i1] = val4 * E[i1-1];
			++i1;//update i1 by 1
		}

		//storing the ai values
		i2=1;
		while(i2<=P)
		{
			A[f1][i2]= Alpha[P][i2];
			++i2;//update i2 by 1
		}
		++f1;//update f1 by 1
	}
	
	//we are finding cepstral cofficients
	double sum1=0;
	int f2=0,m1=0,k1=0;
	double log_input=0.0,temp;
	while(f2<F)
	{
		log_input= R[f2][0]*R[f2][0];
		C[f2][0]=log(log_input);
		m1=1;
		while(m1<=P)
		{
			sum1=0;
			k1=1;
			while(k1<m1)
			{
				temp= k1*C[f2][k1]*A[f2][m1-k1];
				sum1 = sum1 + temp/(m1*1.0);
				++k1;
			}
			C[f2][m1] = A[f2][m1] + sum1;
			++m1;//update m1 by 1
		}
		++f2;//update f2 by 1
	}
	
	//we are applying raised sine window on Ci's
	long double sum2=0.0,sin_input=0.0;
	int f3=0,m3;
	while(f3<F)
	{
		m3=1;
		while(m3<=P)
		{
			sin_input=(PI * m3)/P;//sine input value is stored
			sum2 = (P/2)*sin(sin_input);
			C[f3][m3] = C[f3][m3]* sum2;
			m3++;//update m by 1
		}
		f3++;//update f3 by 1
	}

	
	//we are saving Ci values into a 3D matrix
	int x=0,y=0;
	while(x<F)
	{
		y=0;
		while(y<P)
		{
			allCi[files][x][y+1]=C[x][y+1];
			++y;
		}
		++x;
	}
	++files;
}


/*utility function takes the input and find DCshift and normalization is done on it
  and finding the stable frames and storing the frames into steadyFrames array
  and applying hamming window on all stable frames,then finding Ris*/
void read_file(char *filename)
{
	double mx, PI=3.142857142857; //mx store max value in our sample
	double dcShift, nFactor; //to store Dcshift and normalization value
    double limit = 5000.0;//limit is taken 5000 for normalization
	
    long int enSize; //store size of array x[]
	long int xSize; //store size of array energy[]

	char line1[50];
    FILE *fp = fopen(filename, "r");

    if(fp == NULL)
	{
        printf("unable to open file\n");
		exit(0);
    }

    //to find max value for normalization
    mx = 0;
	long int totalSample1 = 0;
	double temp;
    for(;!feof(fp);)
	{
        fgets(line1, 50, fp);
        if(isalpha(line1[0]) == 0)
		{
            totalSample1++;
			temp = abs(atoi(line1));
			mx=(mx < temp)?temp:mx;
        }
    }
    fclose(fp);
    nFactor = (double)limit/mx;
   
    //doing the Dcshift
    long int sample_count1 = 0;
    FILE *fp1 = fopen(filename, "r");
	char line2[50];
    
    if(fp1 == NULL)
	{
        printf("File not found\n");
        exit(0);
    }
    
	dcShift = 0;
    for(;!feof(fp1);)
	{
        fgets(line2, 50, fp1);
        dcShift = dcShift + atof(line2);
        sample_count1 = sample_count1 + 1;
    }
    dcShift = dcShift / sample_count1;
    fclose(fp1);

	//doing normalization, reading the input values and after normalization storing ,writing in txt file
	char line[50];
	FILE *ip = fopen(filename, "r");
    if(ip == NULL)
	{
	   printf("Error in opening file \n");
	   exit(0);
	}
	xSize = 0;//store 0 initially
	int y;
	double normalizedX,val;
    for(;!feof(ip);)
	{
        fgets(line,50, ip);
        if(isalpha(line[0]) == 0)
		{
            y = atof(line);
			val= (y-dcShift);
            normalizedX = floor(val*nFactor);
            x[xSize] = normalizedX;
			xSize++;
        }
    }
	fclose(ip);

	//now we will find stable frames
	long int totalSample = 0, max_i = 0;
	int n = 0;
	double en = 0, m = 0,val1;
	enSize = 0;
	//to calculate the short term energy
	double energy[80000];

	for(;totalSample < xSize;totalSample++)
	{
		if(n == N)//if condition satisfies, store short term energy to arry energy[]
		{
			en = en / N;//finding average of energy
			if(m < en)// condition for max energy frame or not
			{
			  max_i = enSize;//marking index of highest energy
			  m = en; 
			}

			energy[enSize] = en; //storing short term energy value
			enSize++;
			//reset n and en for next iteration
			n = 0;
			en = 0;
		}
		n++;
		val1 = x[totalSample] * x[totalSample];
		en = en + val1  ;
	}

    /*start_i stores start marker of steady frames
	  and end_i stores end marker of steady frames*/
    long int start_i=0, end_i=0;
	double steadyFrames[F][N];
	if(max_i > 2)//findind start marker
	{
      start_i = (max_i-2)*N;
	}
	else
	{
      start_i=0;
	}

	if(max_i < enSize-3)//findind end marker
	{
      end_i = (max_i+3)*N;
	}
	else
	{
      end_i=enSize*N;
	}

	int f = 0;
	int i = start_i, j=0;
	
	while( i<end_i)
	{
		//storing stable frames into steadyFrames[] array
		steadyFrames[f][j++] = x[i];
		if(j == N)
		{
		  f++;
		  j=0;
		}
		++i;
	}

	//we are applying hamming window to all stable frames
	int i1=0,j1=0;
	double cos_input;
	while(i1<F)
	{
		j1=0;
		while(j1<N)
		{
			cos_input=(2*PI*steadyFrames[i1][j1])/N-1;
			steadyFrames[i1][j1] = steadyFrames[i1][j1] * (0.54-0.46*cos(cos_input));
			++j1;
		}
		++i1;
	}

	//calculating Ris using stable frames
	int f1=0,m1=0,k1=0;
	double x[2]={0.0,0.0};
	while(f1<5)
	{
		m1=0;
		while(m1<=P)
		{
			R[f1][m1] = 0;
			k1=0;
			while(k1<N-m1)
			{
				x[0] = steadyFrames[f1][k1];
				x[1] = steadyFrames[f1][k1+m1];
				R[f1][m1] = R[f1][m1] + x[0]*x[1];
				k1++;
			}
			++m1;
		}
		++f1;
	}

	durbinAlgo();//call durbinAlgo function to evaluate Ais and to do other tasks
}



//utility fucntion to calculate the distance from Ci values of training sets
double tokhuraDistance(FILE *ip)
{
	char line[500];
	int f = 0,size=0;
	double finalDist=0.0;
    double tokhuraWeights[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};//tokhura weights

	for(;!feof(ip) && f<F;)//to read reference data sets and store in array restoreCi
	{
		size= sizeof(line);
		fgets(line,size, ip);
		sscanf(line, "%f %f %f %f %f %f %f %f %f %f %f %f %f", &restoreCi[f][1], &restoreCi[f][2], &restoreCi[f][3], &restoreCi[f][4], &restoreCi[f][5], &restoreCi[f][6], &restoreCi[f][7], &restoreCi[f][8], &restoreCi[f][9], &restoreCi[f][10], &restoreCi[f][11], &restoreCi[f][12]);
		f++;
	}

	int i=0,p=0;
	double d,x,dist,val;
	while(i<F)
	{
		dist = 0;
		p=1;
		while(p<=P)
		{
			d = (C[i][p]- restoreCi[i][p]);
			x=tokhuraWeights[p-1]*d*d;
			dist = dist + x;
			++p;
		}
		val= dist/(P*1.0);
		finalDist = finalDist + val;
		++i;
	}
	double finalval= finalDist/(F*1.0);//evaluating the average
	return finalval; 
}

//utility function for finding the distance
char calculateTokhura()
{
	char predictedVowel;//to store prediction value of vowel
	char filename[30];
	FILE *ip;
	int i1=0,i2=0,p=0,f=0;
	double minDist = 999999;
    double tokhuraWeights[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};//tokhura weights
    char vowels[5] = {'a', 'e', 'i', 'o', 'u'};//to store vowels
	double finalDist,d,dist,distance;
	while(i1<5)
	{
		sprintf(filename, "Reference/reference_ci_%c.txt", vowels[i1]);
		ip = fopen(filename, "r");
		if(ip == NULL)
		{
		 printf("Error in opening file %s\n", filename);
		 exit(0);
		}
		char line[2000];
	    f = 0;
		int size=0;
	   
	    for(;!feof(ip) && f<F;) //store in restoredCi[] after reading the reference file
		{
		size= sizeof(line);
		fgets(line, size, ip);
		sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &restoreCi[f][1], &restoreCi[f][2], &restoreCi[f][3], &restoreCi[f][4], &restoreCi[f][5], &restoreCi[f][6], &restoreCi[f][7], &restoreCi[f][8], &restoreCi[f][9], &restoreCi[f][10], &restoreCi[f][11], &restoreCi[f][12]);
		f++;
	    }
        finalDist = 0;
		i2=0;
	    while(i2<F)
		{
		   dist = 0;
		   p=1;
		   while(p<=P)
		   {
		    d = (C[i2][p]- restoreCi[i2][p]);
			dist = dist + tokhuraWeights[p-1]*d*d;
			++p;
		   }
		 finalDist = finalDist + dist/(P*1.0);
		 ++i2;
	    }
		double val=(F*1.0);
		distance = finalDist/val; // evaluating the average
		tokhuraDist[i1] = distance;//storing the evaluated distance
		
		if(minDist > distance)//we are checking for min distance
		{
			predictedVowel = vowels[i1];//update vowel that we want to returned
			minDist = distance;//update if it is the min distance
		}
	  ++i1;
	}
	return predictedVowel;
}


int _tmain(int argc, _TCHAR* argv[])
{
	//***************************training of test files***************************
    //training using 50 recordings of 10 recording per vowel
	char filenames[5][30] = {"Training/224101044_a_0.txt", "Training/224101044_e_0.txt", "Training/224101044_i_0.txt", "Training/224101044_o_0.txt", "Training/224101044_u_0.txt"};//
	int i=0,j=0;
	while(i<5)//loop over vowels
	{
		j=0;
		while(j<10)//loop over files
		{
			read_file(filenames[i]);
			filenames[i][21]++;//go for next file
			++j;
		}
		++i;
	}

	//***************************generating reference files***************************
	FILE *indi;
	char filename[30];
	int index = 0;
	int vl=0,f=0,p=0,file=0;
	double avgCi[25][P+1];
    char vowels[5] = {'a', 'e', 'i', 'o', 'u'};	//to store vowels
	while(vl<5)//loop over vowel to save file 
	{
		sprintf(filename, "Reference/reference_ci_%c.txt", vowels[vl]);
		indi = fopen(filename, "w");
		f=0;
		while(f<F)//looping over frames
		{
			p=0;
			while(p<P)//looping over p
			{
				double sum = 0;
				//finding average of Cis
				file=vl*10;
				while(file < (vl+1)*10 )
				{
					sum = sum + allCi[file][f][p+1];
					++file;
				}
				sum = sum / 10.0;
				avgCi[index][p+1] = sum;
				fprintf(indi, "%lf ", sum);
				++p;
			}
			index++;
			fprintf(indi, "\n");
			++f;
		}
		fclose(indi);
		vl++;
	}

	//***************************testing of the text files***************************
	char filename1[5][30] = {"Testing/224101044_a_10.txt", "Testing/224101044_e_10.txt", "Testing/224101044_i_10.txt", "Testing/224101044_o_10.txt", "Testing/224101044_u_10.txt"};
	files = 0;//files variable is initialized to 0

	int x=0,file1=0;
    int totalCorrect = 0, individualCorrect= 0;//variable to find accuracy
	double accuracy1,accuracy2;
	while(x<5)
	{
		individualCorrect = 0;
		file1=0;
		while(file1<10)
		{
			read_file(filename1[x]);
			char prediction = calculateTokhura();//calculateTokhura() will give the predicted vowel
			printf("\nfile: %s recognized as = %c\n ", filename1[x],prediction);
			
			//if prediction is accurate then increasing values of individualCorrect & totalCorrect
			individualCorrect = (prediction == vowels[x])?individualCorrect+1:individualCorrect;
			totalCorrect = (prediction == vowels[x])?totalCorrect+1:totalCorrect;

			filename1[x][21]++;//go for next file
			++file1;
		}
		accuracy1= (individualCorrect/10.0)*100;//calculating the accuarcy of a vowel
		printf("**************************** Accuracy of vowel %c = %.3lf %% ****************************\n\n", vowels[x],accuracy1);//printing the accuracy
		++x;
	}
	accuracy2= (totalCorrect/50.0)*100;//calculating the overall accuracy
	printf("************** Overall Accuracy = %.3lf %% **************\n",accuracy2);//printing the overall accuracy
	
	system("PAUSE");
	return 0;
}

