// a111.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


int _tmain(int argc, _TCHAR* argv[])
{  
	int arr[5]={2,4,1,8,3};
	for(int i=0;i<5;i++){
		int flag=-1;
		for(int j=0;j<4;j++){
			if(arr[j]>arr[j+1]){
				int t=arr[j];
				arr[j]=arr[j+1];
				arr[j+1]=t;
				flag++;
			}
		}
		if(flag==-1) break;
	}
	for(int i=0;i<5;i++){
	printf("%d ",arr[i]);
	}
	return 0;
}

