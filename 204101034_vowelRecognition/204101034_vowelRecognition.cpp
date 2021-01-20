// 204101034_vowelRecognition.cpp : Defines the entry point for the console application.


#include "stdafx.h"
#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include<cmath>

using namespace std;

int p=12;
int track=0,count_i=-1;
long double m=0, Ci[50][12], restore_Ci[25][12],test_Ci[5][12];
long double w[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
long double tokh_Dist[25];
int u=0,v=0;

/*this function apply dc shift to the sample values*/
void dc_Shift(vector<long double> *v)
{
	long double sum=0,dc=0;
	int i=0;
	while(i<1000)
	{
		sum += (*v)[i];
		i++;
	}
		
	dc=sum/i;

	
	for(int i=0;i<v->size();i++)
	{
		(*v)[i] -= dc;
	}
		
	if(dc==0)
	{
	//	cout<<"No DC-Shift"<<endl;
	}
	else
	{
	//	cout<<"DC-Shift"<<endl;
	}
}

/*this function normalize the sample value*/
void normalization(vector<long double> *v, long double max)
{
	//cout<<"\n Normaliztion"<<endl;
	long double factor;
	//cout<<"\nMaximum = "<<max<<", Index = "<<track<<endl;
	factor=10000/max;
	//cout<<"\n Factor = "<<factor<<endl;
	for(int i=0;i<v->size();i++)
	{
		((*v)[i])*=factor;
	}


}

/*This function apply Hamming window in a single frame of size 320 samples*/
void hamming_Window(vector<long double> *v)
{
	//cout<<"\n Hamming Window ";
	long double temp=0;
	int N=v->size();

	//cout<<N<<endl;

	for(int i=0;i<v->size();i++)
	{
		temp = 0.54 - (0.46 * cos(2*3.14*i/(N-1)));
		((*v)[i]) *= temp;
	}

	
}

/*This function apply the Raised Sine Window in Ci of each frame*/
void raised_SineWindow(long double *C)
{
	////cout<<"\n after Raised Sine Window"<<endl;
	long double sum=0;
	for(int i=1;i<=p;i++)
	{
		sum= (p/2)*	sin((3.14*C[i])/p);
		C[i]=1 + sum ;
		////cout<<"C["<<i<<"]="<<C[i]<<endl;	
	}
}


/*This function find the maximum sample value and return it*/
long double find_Max(long double x)
{
	count_i++;
	if(x<0)
	{
		x*=(-1);
	}
	if(m>x)
	{
		m=m;
	}
	else
	{
		m=x;
		track=count_i;
	}
	return m;
}

/*This Function Create Average Ci for each vowel for 5 Frame each*/
void avg_Ci(string s)
{
	ifstream in;
	ofstream out;
	out.open(s.c_str());
	vector<vector<long double> > v;
	long double Avg[5][12];
	long double x;
	int i,j;

	for(int i=0;i<5;i++)
	{
		for(int j=0;j<12;j++)
			Avg[i][j]=0;
	}
	for(int i=0,m=0;i<50,m<5;)
	{
		
		for(int j=0,n=0;j<12,n<12;j++,n++)
		{
			Avg[m][n]+=Ci[i][j];
			
			if(i==45+m)
			Avg[m][n]/=10;
		}
		
		i++;
		m++;
		if(m==5 && i!=50)
			m=0;
		
	}
	 
	 for( i=0;i<5;i++)
	{
		for( j=0;j<12;j++)
		{	
			out<<Avg[i][j];
			out<<" ";
		}
		out<<endl;

	}
	out.close();
	
}


/*This Function Help to hold Ci values in a 2d Array*/
void store(long double * C)
{
	for(v=0;v<12;v++)
		Ci[u][v]=C[v+1];
		
	u++;
}


/*This function calulate the cepstral coeff Ci's*/
void cepstral_Coeff(long double *R,long double *a)
{
	long double C[13];
	long double sum=0;
	
	C[0]=log(R[0]);
	
	for(int m=1;m<=p;m++)
	{
		sum=0;
		for(int k=1;k<=m-1;k++)
		{
			sum += (k*C[k]*a[m-k])/m;
		}
		C[m]=a[m]+sum;
		
	
	}
	raised_SineWindow(C);
	
	store(C);
}

/*This Function Apply Durbin Algorithm And Find The value of ai's */
void durbin(long double *R)
{
	long double Alpha[13][13],E[13],K[13],a[13];
	long double sum=0;

	E[0]=R[0];
	
	for(int i=1;i<=p;i++)
	{
		sum=0;
		for(int j=1;j<=i-1;j++)
		{
			sum += Alpha[i-1][j]*R[i-j];	
		}
			
		K[i]=(R[i]-sum)/E[i-1];
				
		Alpha[i][i]=K[i];
		
		for(int j=1;j<=i-1;j++)
		{
			Alpha[i][j]=Alpha[i-1][j] - K[i]*Alpha[i-1][i-j];
		}
		
		E[i]=(1-pow(K[i],2))*E[i-1];
		
		
	}
	for(int i=0;i<=p;i++)
	{
		a[i]=Alpha[p][i];
	}

	cepstral_Coeff(R,a);
	
}

/*This Function Evaluate the Values of Ri's */
void autocorrelation(vector<long double> *v)
{
	long double R[13];
	
	for(int i=0;i<=12;i++)
	{
		R[i]=0;
		for(int j=0;j<320-i;j++)
		{
			R[i]+= ((*v)[j])*((*v)[j+i]);
		}
	}
	
	durbin(R);
}


/*This function make the frame of size 320 sample*/
void framing(vector<long double> *v, int track)
{
	vector<long double> temp;
	int start,end,count=0,j=0;
	if(track >4000)
	{
		start=track-800;
		end=track+799;
	}
	else
	{
		start=track;
		end=track+1599;
	}
	
	for(int i=start;i<=end;i++)
	{

		temp.push_back((*v)[i]);
		count++;
		if(count==320)
		{
			j++;
			
			hamming_Window(&temp);
			
			autocorrelation(&temp);

			count=0;
			
			temp.clear();//this clear the value of vector temp	
		}
		
		
	}
	
}

/*This function Read the text File and call fuction which evaluate further process*/
void read_File(string s)
{
	ifstream in1;
	ofstream out1;
	int i=0;
	in1.open(s.c_str());
	long double max,input;
	vector<long double> v1;
	if(in1.is_open())
	{
		while(!in1.eof())
		{
			in1>>input;
			max=find_Max(input);
			v1.push_back(input);
			
		}
	}
	else
	{
		cout<<s<<" file not open"<<endl;
		exit(1);
	}
	
	dc_Shift(&v1);
	
	normalization(&v1,max);
	
	framing(&v1,track);	

}

/*This function create the Reference File*/
void make_RefFile()
{
	string training_File[5]={"204101034_a_0.txt","204101034_e_0.txt","204101034_i_0.txt","204101034_o_0.txt","204101034_u_0.txt"};
	
	for(int i=0;i<5;i++)
	{
		for(int j=0;j<10;j++)
		{
			//cout<<"file("<<i<<","<<j<<")"<<endl;
			
			m=0;count_i=-1;
			
			//cout<<training_File[i]<<endl;
			
			read_File(training_File[i]);
			
			training_File[i][12]++;
			
		}
		if(i==0)
		{
			avg_Ci("Avg_ref_a1.txt");
			cout<<"Avg_ref_a1.txt"<<endl;
		}
		if(i==1)
		{
			avg_Ci("Avg_ref_e1.txt");
			cout<<"Avg_ref_e1.txt"<<endl;
		}
		if(i==2)
		{	
			avg_Ci("Avg_ref_i1.txt");
			cout<<"Avg_ref_i1.txt"<<endl;
		}
		if(i==3)
		{
			avg_Ci("Avg_ref_o1.txt");
			cout<<"Avg_ref_o1.txt"<<endl;
		}
		if(i==4)
		{
			avg_Ci("Avg_ref_u1.txt");
			cout<<"Avg_ref_u1.txt"<<endl;
		}
			
		u=0;
		v=0;
	}
}

/*This Function restore the Reference file for testing process*/
void restore_RefFile(string s)
{
	ifstream in1(s.c_str());
	long double x;
	if(in1.is_open())
	{
		do
		{
			for(v=0;v<12;v++)
			{
				in1>>x;
				restore_Ci[u][v]=x;
			}

			u++;
		}while(u%5!=0);
	}
	else
	{
		cout<<"file not open"<<endl;
	}
	
}

/*This Algorithm find the Tokhura Distance and Tell the which vowel is recognized by finding the minimum distance*/
void tokhuraAlgo()
{
	//cout<<"Tokhura\n";
	long double min=0;
	int index;
	int m=0;
	
	for(int i=0;i<25,m<5;)
	{
		tokh_Dist[i]=0;
		for(int j=0;j<12;j++)
		{
			tokh_Dist[i]=tokh_Dist[i]+ w[j]*pow((Ci[m][j]-restore_Ci[i][j]),2);
		}
		
		i++;
		m++;
		
		if(m==5 && i!=25)
		{
			m=0;
		}	
			
	}

	min= tokh_Dist[0];
	index=0;

	for(int i=0;i<25;i++)
	{
		if(min>tokh_Dist[i])
		{
			min = tokh_Dist[i];
			index = i;
		}	
	}
	
	
	if((index >= 0) && (index <= 4))
	cout<<"Vowel a"<<endl;
	
	if((index >= 5) && (index <= 9))
	cout<<"Vowel e"<<endl;
	
	if((index >= 10) && (index <= 14))
	cout<<"Vowel i"<<endl;
	
	if((index >= 15) && (index <= 19))
	cout<<"Vowel o"<<endl;
	
	if((index >= 20) && (index <= 24))
	cout<<"Vowel u"<<endl;
	
}

/*This function Read the Test file and Process the data and call functions which Evaluate the Ci's for each frame*/
void read_TestFile(string s)
{
	ifstream in1;
	ofstream out1;
	int i=0;
	long double max,input;
	vector<long double> v1;


	in1.open(s.c_str());
	
	if(in1.is_open())
	{
		while(!in1.eof())
		{
			in1>>input;
			
			max=find_Max(input);
			
			v1.push_back(input);
		}
	}
	else
	{
		cout<<s<<" file not open"<<endl;
		exit(1);
	}

	dc_Shift(&v1);
	
	normalization(&v1,max);
	
	framing(&v1,track);	
	
}
 

/*This Fuction called two function 1st evaluate the Ci's of Test File and 2nd Calculate the tokhura distance*/
void cal_TokhuraDistance(string s)
{
	read_TestFile(s);
	
	tokhuraAlgo();
}

/*This Function is used to testing from the text file*/
void testing_Text()
{
	ifstream in1,in2,in3,in4,in5,in6;
	ofstream out;
	u=0,v=0;
		
	restore_RefFile("Avg_ref_a1.txt");
	restore_RefFile("Avg_ref_e1.txt");
	restore_RefFile("Avg_ref_i1.txt");
	restore_RefFile("Avg_ref_o1.txt");
	restore_RefFile("Avg_ref_u1.txt");
		
	u=0;

	string testing_File[5]={"204101034_a_10.txt","204101034_e_10.txt","204101034_i_10.txt","204101034_o_10.txt","204101034_u_10.txt"};
	
	for(int i=0;i<5;i++)
	{
		for(int j=0;j<10;j++)
		{
			//cout<<"file("<<i<<","<<j<<")"<<"  =>";
			m=0;count_i=-1;
			
			cout<<" File : "<<testing_File[i]<<" => Recognized ";
			cal_TokhuraDistance(testing_File[i]);
			u=0;
			testing_File[i][13]++;
			
		}
		cout<<endl;
	}

}

/*This Function is for live testing*/
void testing_Mic()
{
	ifstream in1,in2,in3,in4,in5,in6;
	ofstream out;
	u=0,v=0;
		
	restore_RefFile("Avg_ref_a1.txt");
	restore_RefFile("Avg_ref_e1.txt");
	restore_RefFile("Avg_ref_i1.txt");
	restore_RefFile("Avg_ref_o1.txt");
	restore_RefFile("Avg_ref_u1.txt");
	
	system("Recording_Module.exe 2 input_file.wav input_file.txt");
	
	cal_TokhuraDistance("input_file.txt");
	u=0;
}

int _tmain(int argc, _TCHAR* argv[])
{
	char ch;
	int z;
	cout<<"Making Reference File: \n"<<endl;
	make_RefFile();

	cout<<"\n\nTesting"<<endl;
	do
	{
		cout<<"1.Test by input Text file\n2.Test by microphone\n";
		cin>>z;
		switch(z)
		{
			case 1:testing_Text();
				break;
				
			case 2:testing_Mic();
				break;
				
				
			default:cout<<"Invalid Input";
				break;
		}
		cout<<"do want to continue y/n";
		cin>>ch;
	}while(ch=='y' || ch == 'Y');
	
	return 0;
}
