#include "mdls_pgm.h"

int read_pnm_header (FILE * infile, int *colsp, int *rowsp, int *maxsp,
		 int *typep)
{
  int rows=-1, cols=-1, maxs=-1, type=-1;
  int res;
  char line[LINELEN + 1];
 
  while (fgets (line, LINELEN, infile) != NULL
	 && (line[0] == '#' || line[0] == 0 || line[0] == '\n'));

  res = sscanf (line, "P%d %d %d %d", &type,&cols,&rows,&maxs);
  if ((res < 1) || (type != 5 && type != 6 && type != 7))
    {
      fprintf (stderr, "bad header: type\n");
      return -1;
    }
  if (res < 4) { /* read header, need cols rows and maxs */
    while (fgets (line, LINELEN, infile) != NULL
  	 && (line[0] == '#' || line[0] == 0 || line[0] == '\n'));
    res = sscanf (line, "%d %d %d", &cols, &rows,&maxs);
    if (res < 2)
    {
      fprintf (stderr, "bad header: cols rows\n");
      return -1;
    }
    if (res < 3)  { /* read header, cols, rows; need maxs */
      while (fgets (line, LINELEN, infile) != NULL
  	 && (line[0] == '#' || line[0] == 0 || line[0] == '\n'));
      if (sscanf (line, "%d", &maxs) != 1)
      {
        fprintf (stderr, "bad header: maxs\n");
        return -1;
      }
    }
  } 
  *colsp = cols;
  *rowsp = rows;
  *maxsp = maxs;
  *typep = type;
  return 0;
}
