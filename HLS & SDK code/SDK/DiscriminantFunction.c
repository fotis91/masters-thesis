/******************************************************************************
*
* Copyright (C) 2009 - 2014 Xilinx, Inc.  All rights reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* Use of the Software is limited solely to applications:
* (a) running on a Xilinx device, or
* (b) that interact with a Xilinx device through a bus or interconnect.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
* XILINX  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF
* OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
* Except as contained in this notice, the name of the Xilinx shall not be used
* in advertising or otherwise to promote the sale, use or other dealings in
* this Software without prior written authorization from Xilinx.
*
******************************************************************************/

/*
 * Author: Fotios Kostarelos
 *
 * This application is the driver code in C that controls the HlS
 * hardware accelerators, thus implementing the discriminant  
 * function for TBI detection.
 *
 *   
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "platform.h"
#include "xil_printf.h"
#include "xparameters.h"
#include "xaxidma.h"
#include "xamplitudeasymmetry.h"
#include "xcoypd.h"
#include "xrelativepower.h"



void initializations(XAxiDma *axiDma1, XAxiDma *axiDma2, XAxiDma *axiDma3,XAxiDma *axiDma4,XAxiDma *axiDma5,
		XAmplitudeasymmetry *XAmplitudeasymmetryInitPtr, XCoypd *XCoypdInitPtr,XRelativepower *XRelativepowerInitPtr);
int init_dma(XAxiDma *axiDma,int dmaid);
int amplitude_asymmetry_init(XAmplitudeasymmetry *XAmplitudeasymmetryInitPtr);
int coypd_init(XCoypd *XCoypdInitPtr);
int relativepower_init(XRelativepower *XRelativepowerInitPtr);

void Double_Channel_transaction(float * A, float * B, XAxiDma *axiDma1,XAxiDma *axiDma2);
void Single_Channel_transaction(float * A, XAxiDma *axiDma1);


#define transfersize 8192

float P3[transfersize] = {);
 
float P4[transfersize] = {};
 
float O1[transfersize] = {};

float O2[transfersize] = {};
 
float T8[transfersize] = {};

float P7[transfersize] = {};

float P8[transfersize] = {};
 
float F3[transfersize] = {};
 
float F7[transfersize] = {};

float Fp1[transfersize] = {};

float T7[transfersize] = {};
 
float C3[transfersize] = {};

float Fp2[transfersize] = {};

float F4[transfersize] = {};
 
float F8[transfersize] = {};

int main()
{
	//variable declarations
	XAxiDma axiDma1,axiDma2,axiDma3,axiDma4,axiDma5;
	XAmplitudeasymmetry amplitudeasymmetry;
	XCoypd coypd;
	XRelativepower relativepower;
	int result=0,result1=0;
	
	float Alpha_Amplitude_Asymmetry_between_F3_O1, Alpha_Amplitude_Asymmetry_between_F4_O2, Alpha_Amplitude_Asymmetry_between_F4_P8,
	Alpha_Amplitude_Asymmetry_between_F8_P8, Alpha_Amplitude_Asymmetry_between_O1_F7, Alpha_Relative_Power_for_O1,
	Alpha_Relative_Power_for_O2, Alpha_Relative_Power_for_P3, Alpha_Relative_Power_for_P4, Alpha_Relative_Power_for_P7,
	Alpha_Relative_Power_for_P8, Alpha_Relative_Power_for_T8, Beta_Amplitude_Asymmetry_between_F4_O2, Beta_Amplitude_Asymmetry_between_F4_P8,
	Beta_Amplitude_Asymmetry_between_F8_P8, Beta_Coherence_between_C3_P3, Beta_Coherence_between_T7_P7, Beta_Phase_Difference_between_F3_F4,
	Beta_Phase_Difference_between_Fp2_F4, Theta_Coherence_between_Fp1_F3,Beta_Amplitude_Asymmetry_between_F3_O1,
	Beta_Amplitude_Asymmetry_between_O1_F7;

	//components initializations
	initializations(&axiDma1,&axiDma2,&axiDma3,&axiDma4,&axiDma5,&amplitudeasymmetry,&coypd,&relativepower);

	//calculation of the 20 variables for the DF as described in the thesis
	////////////COHERENCE FOR CHANNELS Fp1-F3
	XCoypd_Set_sel1(&coypd,1);
	XCoypd_Set_sel2(&coypd,1);
	XCoypd_Start(&coypd);
	Double_Channel_transaction(Fp1,F3,&axiDma2,&axiDma3);
	result = XCoypd_Get_output_r(&coypd);
	Theta_Coherence_between_Fp1_F3=*((float*)&result);

    ////////////COHERENCE FOR CHANNELS T7-P7
	XCoypd_Set_sel1(&coypd,1);
	XCoypd_Set_sel2(&coypd,2);
	XCoypd_Start(&coypd);
	Double_Channel_transaction(T7,P7,&axiDma2,&axiDma3);
	result = XCoypd_Get_output_r(&coypd);
	Beta_Coherence_between_T7_P7=*((float*)&result);

	////////////COHERENCE FOR CHANNELS C3-P3
	XCoypd_Set_sel1(&coypd,1);
	XCoypd_Set_sel2(&coypd,2);
	XCoypd_Start(&coypd);
	Double_Channel_transaction(C3,P3,&axiDma2,&axiDma3);
	result = XCoypd_Get_output_r(&coypd);
	Beta_Coherence_between_C3_P3=*((float*)&result);

	////////////PHASE DIFFERNCE FOR CHANNELS Fp2-F4
	XCoypd_Set_sel1(&coypd,2);
	XCoypd_Set_sel2(&coypd,0);
	XCoypd_Start(&coypd);
	Double_Channel_transaction(Fp2,F4,&axiDma2,&axiDma3);
	result = XCoypd_Get_output_r(&coypd);
	Beta_Phase_Difference_between_Fp2_F4=*((float*)&result);

	////////////PHASE DIFFERNCE FOR CHANNELS F3-F4
	XCoypd_Set_sel1(&coypd,2);
	XCoypd_Set_sel2(&coypd,0);
	XCoypd_Start(&coypd);
	Double_Channel_transaction(F3,F4,&axiDma2,&axiDma3);
	result = XCoypd_Get_output_r(&coypd);
	Beta_Phase_Difference_between_F3_F4=*((float*)&result);


	////////////AMPLITUDE ASYMMETRY FOR CHANNELS F4-P8
	XAmplitudeasymmetry_Start(&amplitudeasymmetry);
	Double_Channel_transaction(F4,P8,&axiDma4,&axiDma5);
	result = XAmplitudeasymmetry_Get_AA1(&amplitudeasymmetry);
	Alpha_Amplitude_Asymmetry_between_F4_P8=*((float*)&result);
	result1 = XAmplitudeasymmetry_Get_AA2(&amplitudeasymmetry);
	Beta_Amplitude_Asymmetry_between_F4_P8=*((float*)&result1);

	////////////AMPLITUDE ASYMMETRY FOR CHANNELS F8-P8
	XAmplitudeasymmetry_Start(&amplitudeasymmetry);
	Double_Channel_transaction(F8,P8,&axiDma4,&axiDma5);
	result = XAmplitudeasymmetry_Get_AA1(&amplitudeasymmetry);
	Alpha_Amplitude_Asymmetry_between_F8_P8=*((float*)&result);
	result1 = XAmplitudeasymmetry_Get_AA2(&amplitudeasymmetry);
	Beta_Amplitude_Asymmetry_between_F8_P8=*((float*)&result1);

	////////////AMPLITUDE ASYMMETRY FOR CHANNELS F4-O2
	XAmplitudeasymmetry_Start(&amplitudeasymmetry);
	Double_Channel_transaction(F4,O2,&axiDma4,&axiDma5);
	result = XAmplitudeasymmetry_Get_AA1(&amplitudeasymmetry);
	Alpha_Amplitude_Asymmetry_between_F4_O2=*((float*)&result);
	result1 = XAmplitudeasymmetry_Get_AA2(&amplitudeasymmetry);
	Beta_Amplitude_Asymmetry_between_F4_O2=*((float*)&result1);

	////////////AMPLITUDE ASYMMETRY FOR CHANNELS F3-O1
	XAmplitudeasymmetry_Start(&amplitudeasymmetry);
	Double_Channel_transaction(F3,O1,&axiDma4,&axiDma5);
	result = XAmplitudeasymmetry_Get_AA1(&amplitudeasymmetry);
	Alpha_Amplitude_Asymmetry_between_F3_O1=*((float*)&result);
	result1 = XAmplitudeasymmetry_Get_AA2(&amplitudeasymmetry);
	Beta_Amplitude_Asymmetry_between_F3_O1=*((float*)&result1);

	////////////AMPLITUDE ASYMMETRY FOR CHANNELS O1-F7
	XAmplitudeasymmetry_Start(&amplitudeasymmetry);
	Double_Channel_transaction(O1,F7,&axiDma4,&axiDma5);
	result = XAmplitudeasymmetry_Get_AA1(&amplitudeasymmetry);
	Alpha_Amplitude_Asymmetry_between_O1_F7=*((float*)&result);
	result1 = XAmplitudeasymmetry_Get_AA2(&amplitudeasymmetry);
	Beta_Amplitude_Asymmetry_between_O1_F7=*((float*)&result1);


	//Relative power channel P3
	XRelativepower_Start(&relativepower);
	Single_Channel_transaction(P3,&axiDma1);
	result = XRelativepower_Get_RP(&relativepower);
	Alpha_Relative_Power_for_P3=*((float*)&result);

	//Relative power channel P4
	XRelativepower_Start(&relativepower);
	Single_Channel_transaction(P4,&axiDma1);
	result = XRelativepower_Get_RP(&relativepower);
	Alpha_Relative_Power_for_P4=*((float*)&result);

	//Relative power channel O1
	XRelativepower_Start(&relativepower);
	Single_Channel_transaction(O1,&axiDma1);
	result = XRelativepower_Get_RP(&relativepower);
	Alpha_Relative_Power_for_O1=*((float*)&result);

	//Relative power channel O2
	XRelativepower_Start(&relativepower);
	Single_Channel_transaction(O2,&axiDma1);
	result = XRelativepower_Get_RP(&relativepower);
	Alpha_Relative_Power_for_O2 = *((float*)&result);

	//Relative power channel T8
	XRelativepower_Start(&relativepower);
	Single_Channel_transaction(T8,&axiDma1);
	result = XRelativepower_Get_RP(&relativepower);
	Alpha_Relative_Power_for_T8 = *((float*)&result);

	//Relative power channel P7
	XRelativepower_Start(&relativepower);
	Single_Channel_transaction(P7,&axiDma1);
	result = XRelativepower_Get_RP(&relativepower);
	Alpha_Relative_Power_for_P7 = *((float*)&result);

	//Relative power channel P8
	XRelativepower_Start(&relativepower);
	Single_Channel_transaction(P8,&axiDma1);
	result = XRelativepower_Get_RP(&relativepower);
	Alpha_Relative_Power_for_P8 = *((float*)&result);


	return 0;
}

//This function performs a signle transaction for one DMA
void Single_Channel_transaction(float * A, XAxiDma *axiDma1){

	int status;

	Xil_DCacheFlushRange((unsigned)A, transfersize * sizeof(float));

	status = XAxiDma_SimpleTransfer(axiDma1, (u32)A, transfersize * sizeof(float), XAXIDMA_DMA_TO_DEVICE);

	do {
		status = XAxiDma_Busy(axiDma1, XAXIDMA_DMA_TO_DEVICE);
	}while(status);

}

//This function performs a signle transaction for two DMAs
void Double_Channel_transaction(float * A, float * B,XAxiDma *axiDma1,XAxiDma *axiDma2){

	int status;

	Xil_DCacheFlushRange((unsigned)A, transfersize * sizeof(float));
	Xil_DCacheFlushRange((unsigned)B, transfersize * sizeof(float));

	status = XAxiDma_SimpleTransfer(axiDma1, (u32)A, transfersize * sizeof(float), XAXIDMA_DMA_TO_DEVICE);
	status = XAxiDma_SimpleTransfer(axiDma2, (u32)B, transfersize * sizeof(float), XAXIDMA_DMA_TO_DEVICE);

	do {
		status = XAxiDma_Busy(axiDma2, XAXIDMA_DMA_TO_DEVICE);
	}while(status);

}

//This function initializes the DMAs
int init_dma(XAxiDma *axiDmaPtr,int dmaid){
   XAxiDma_Config *CfgPtr;
   int status;
   // Get pointer to DMA configuration
   if(dmaid==1) CfgPtr = XAxiDma_LookupConfig(XPAR_AXIDMA_0_DEVICE_ID);
   else if(dmaid==2) CfgPtr = XAxiDma_LookupConfig(XPAR_AXIDMA_1_DEVICE_ID);
   else if(dmaid==3) CfgPtr = XAxiDma_LookupConfig(XPAR_AXIDMA_2_DEVICE_ID);
   else if(dmaid==4) CfgPtr = XAxiDma_LookupConfig(XPAR_AXIDMA_3_DEVICE_ID);
   else CfgPtr = XAxiDma_LookupConfig(XPAR_AXIDMA_4_DEVICE_ID);


   if(!CfgPtr){
      print("Error looking for AXI DMA config\n\r");
      return XST_FAILURE;
   }
   // Initialize the DMA handle
   status = XAxiDma_CfgInitialize(axiDmaPtr,CfgPtr);
   if(status != XST_SUCCESS){
      print("Error initializing DMA\n\r");
      return XST_FAILURE;
   }
   //check for scatter gather mode - this example must have simple mode only
   if(XAxiDma_HasSg(axiDmaPtr)){
      print("Error DMA configured in SG mode\n\r");
      return XST_FAILURE;
   }
   //disable the interrupts
   XAxiDma_IntrDisable(axiDmaPtr, XAXIDMA_IRQ_ALL_MASK,XAXIDMA_DEVICE_TO_DMA);
   XAxiDma_IntrDisable(axiDmaPtr, XAXIDMA_IRQ_ALL_MASK,XAXIDMA_DMA_TO_DEVICE);

   return XST_SUCCESS;
}

//This function initializes the IP for the relativer power calculation
int relativepower_init(XRelativepower *XRelativepowerInitPtr)
{
   XRelativepower *cfgPtr;
   int status;

   cfgPtr = XRelativepower_LookupConfig(XPAR_RELATIVEPOWER_0_DEVICE_ID);
   if (!cfgPtr) {
      print("ERROR: Lookup of acclerator configuration failed.\n\r");
      return XST_FAILURE;
   }
   status = XRelativepower_CfgInitialize(XRelativepowerInitPtr, cfgPtr);
   if (status != XST_SUCCESS) {
      print("ERROR: Could not initialize accelerator.\n\r");
      return XST_FAILURE;
   }
   return status;
}

//This function initializes the IP for the amplitude asymmetry calculation
int amplitude_asymmetry_init(XAmplitudeasymmetry *XAmplitudeasymmetryInitPtr)
{
   XAmplitudeasymmetry *cfgPtr;
   int status;

   cfgPtr = XAmplitudeasymmetry_LookupConfig(XPAR_XAMPLITUDEASYMMETRY_0_DEVICE_ID);
   if (!cfgPtr) {
      print("ERROR: Lookup of acclerator configuration failed.\n\r");
      return XST_FAILURE;
   }
   status = XAmplitudeasymmetry_CfgInitialize(XAmplitudeasymmetryInitPtr, cfgPtr);
   if (status != XST_SUCCESS) {
      print("ERROR: Could not initialize accelerator.\n\r");
      return XST_FAILURE;
   }
   return status;
}

//This function initializes the IP for the coherence and phase difference calculation
int coypd_init(XCoypd *XcoypdInitPtr)
{
   XCoypd *cfgPtr;
   int status;

   cfgPtr = XCoypd_LookupConfig(XPAR_COYPD_0_DEVICE_ID);
   if (!cfgPtr) {
      print("ERROR: Lookup of acclerator configuration failed.\n\r");
      return XST_FAILURE;
   }
   status = XCoypd_CfgInitialize(XcoypdInitPtr, cfgPtr);
   if (status != XST_SUCCESS) {
      print("ERROR: Could not initialize accelerator.\n\r");
      return XST_FAILURE;
   }
   return status;
}

//This function handles the initializations
void initializations(XAxiDma *axiDma1, XAxiDma *axiDma2, XAxiDma *axiDma3,XAxiDma *axiDma4,XAxiDma *axiDma5,
		XAmplitudeasymmetry *amplitudeasymmetry, XCoypd *coypd,XRelativepower *relativepower){
	int status;

	status  = init_dma(axiDma1,1);
		if (status != XST_SUCCESS) {
		 print("dma1 init peripheral setup failed\n\r");
		 exit(-1);
	}

	 status  = init_dma(axiDma2,2);
		if (status != XST_SUCCESS) {
		 print("dma2 init peripheral setup failed\n\r");
		 exit(-1);
	}

	status = init_dma(axiDma3,3);
	  if (status != XST_SUCCESS) {
		print("dma3 init peripheral setup failed\n\r");
		exit(-1);
	}

    status = init_dma(axiDma4,4);
	   if (status != XST_SUCCESS) {
		print("dma4 init peripheral setup failed\n\r");
		exit(-1);
	}

	 status  = init_dma(axiDma5,5);
	  if (status != XST_SUCCESS) {
		print("dma5 init peripheral setup failed\n\r");
		exit(-1);
	}

    status  = relativepower_init(relativepower);
	   if (status != XST_SUCCESS){
		print("relativepower init peripheral setup failed\n\r");
		exit(-1);
	}

   status  = amplitude_asymmetry_init(amplitudeasymmetry);
	   if (status != XST_SUCCESS){
		print("amplitud easymmetry init peripheral setup failed\n\r");
		exit(-1);
	}

   status = coypd_init(coypd);
	   if (status != XST_SUCCESS){
	  print("coherence init peripheral setup failed\n\r");
	   exit(-1);
	}


}


