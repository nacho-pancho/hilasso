/*
 * Matlab C MEX function for simultaneously writing to a log file and printing to matlab console
 *
 * Author  : Ignacio Ramirez  <nacho@fing.edu.uy>
 *
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "mex.h"

#define DEBUG 0

#define MAXLINE 4096
#define TMPLEN 1024

static FILE* LogFileHandle = NULL;

/* 
 * perform basic C-like unescaping since Matlab sends strings 'as is' with the '\t' as two chars
 *
 */
void unescape(char* s, char *d) {
  int i = 0, j = 0;
  while (s[i] != 0) {
    if (s[i] != '\\') {
      d[j++]=s[i++];
    } else {
      switch (s[++i]) {
	case 'n':
	  d[j++]='\n'; i++;
	  break;
	case 't':
	  d[j++]='\t'; i++;
	  break;
	case 'b':
	  d[j++]='\b'; i++;
	  break;
	case '\\':
	  d[j++]='\\'; i++;
	  break;
	case 'r':
	  d[j++]='\r'; i++;
	  break;
	default:
	  d[j++]='\\';
      }
    }
  }
  d[j]='\0';
}

/**
 * 
 *
 * inputs:
 *
 * 0 S ....... printf format string or the reserved strings 'open' or 'close'
 * 1 ......... file name for 'open', ignored for 'close' or first argument to printf format string
 * 2,3... .... further arguments to the format string
 *
 * outputs:
 * n ......... number of chars written
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char formatString[MAXLINE+1];
  char fname[256];
  char *nextField,*prevField;
  int parIndex = 0; /* index to next matlab prhs */
  int i = 0; /* auxiliary index */
  char tmpout[TMPLEN+1];
  char tmp2[MAXLINE+1];
  char fmtToken[32]; /* printf format token */
  char c;
  int L;
  int n;
  double tmpnum;
  double *tmp;
  int tmpRows,tmpCols;
  
  if (nrhs < 1) {
    mexErrMsgTxt("Not enough parameters.");
    return;
  }
  if (mxGetClassID(prhs[0]) != mxCHAR_CLASS) {
    mexErrMsgTxt("First parameter must be a string.");
    return;
  }
  L = mxGetM(prhs[0])*mxGetN(prhs[0])*sizeof(mxChar)+1;
  mxGetString(prhs[0],formatString,L);

  /*
   * open/close
   */
  if (strcmp(formatString,"open") == 0) {
    if (mxGetClassID(prhs[1]) != mxCHAR_CLASS) {
      mexErrMsgTxt("Second parameter must be file name.");
      return;
    }
    L = mxGetM(prhs[1])*mxGetN(prhs[1])*sizeof(mxChar)+1;
    mxGetString(prhs[1],fname,L);

    LogFileHandle = fopen(fname,"w");
    if (LogFileHandle == NULL) {
      mexWarnMsgTxt("Error opening log file. No log will be saved!");
    }
    return;
  } else if (strcmp(formatString,"close") == 0) {
    if (LogFileHandle != NULL) {
      fclose(LogFileHandle);
      LogFileHandle = NULL;
    } else {
      mexWarnMsgTxt("Log file already closed.");      
    }
    return;
  }
  /*
   * write mode: go over format string and fetch parameters as they appear
   */
  if (LogFileHandle == NULL) {
      mexWarnMsgTxt("Log file not open. Nothing will be written.");
  }
  prevField = formatString;
  parIndex = 1;
  do {
    nextField = strchr(prevField,'%');
    if (nextField != NULL) {
      /*
       * fetch next field
       */
      fmtToken[0]='%';
      i = 1;
      while ((c=nextField[i]) != '\0') {
	if (i >= 31) { 
	  break;
	}
        fmtToken[i++]=c;
        if (strchr("dulsfceEgGxXopni%",c) != NULL) { /* end of token */
          break;
        }
      }
      fmtToken[i]=0;
      /* printf("%d:%s\n",parIndex,fmtToken); */
      /* char type */
      if ((c != '%') && (parIndex >= nrhs)) {
        mexWarnMsgTxt("Not enough parameters. output will be truncated.");
	nextField = NULL;
	parIndex++;
        continue;
      }
      switch (c) {
        case 's':
          L = mxGetM(prhs[parIndex])*mxGetN(prhs[parIndex])*sizeof(mxChar)+1;
          mxGetString(prhs[parIndex],tmpout,L);
          break;
        case 'c':
          L = 2;
          mxGetString(prhs[parIndex],tmpout,L);
          break;
	case 'f': case 'g': case 'G': case 'e' : case 'E':
	  tmpnum = (double) *mxGetPr(prhs[parIndex]);
	  if (isnan(tmpnum)) {
	    snprintf(tmpout,TMPLEN,"NaN");
	  }  else if (isinf(tmpnum)) {
	    snprintf(tmpout,TMPLEN,"Inf");
	  } else {
	    snprintf(tmpout,TMPLEN,fmtToken,tmpnum);          
	  }
	  break;
	case 'd': case 'i': case 'l' :
          snprintf(tmpout,TMPLEN,fmtToken,(int) *mxGetPr(prhs[parIndex]));          
	  break;
	case 'u': 
          snprintf(tmpout,TMPLEN,fmtToken,(unsigned int) *mxGetPr(prhs[parIndex]));          
	  break;
	case '%':
          tmpout[0]='%'; tmpout[1]='\0';
	  break;
      }
      /*
       * write actual output
       */
      *nextField = '\0';
      unescape(prevField,tmp2);
      mexPrintf(tmp2);
      mexPrintf(&tmpout[0]);
      if (LogFileHandle != NULL) {
        n += fputs(tmp2, LogFileHandle);
        n += fputs(&tmpout[0], LogFileHandle);
      }
      /* 
       * advance to next argument
       */
      prevField = &nextField[i];
    } else {
      unescape(prevField,tmp2);
      mexPrintf(tmp2);
      if (LogFileHandle != NULL)
	fputs(tmp2, LogFileHandle);
      else
	mexWarnMsgTxt("Log file not open. Nothing will be written.");	
    }
    parIndex++;
  } while (nextField  != NULL);
  fflush(NULL);
}

