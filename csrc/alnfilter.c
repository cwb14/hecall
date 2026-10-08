/* alnfilter: copy the alignment records of a .1aln whose byte in a keep mask is non-zero.
 *
 *   alnfilter <in.1aln> <keep.mask> <out.1aln>
 *
 * keep.mask holds one byte per alignment, in file order. Header lines (sources, skeleton,
 * trace spacing) are copied unchanged, so the output is a valid .1aln for ALNchain / ALNtoPAF.
 * Built against the FASTGA sources (Myers; ONElib by Durbin & Myers); see build.sh.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "GDB.h"
#include "align.h"
#include "alncode.h"

int main(int argc, char *argv[])
{ int64 novl, i, kept;
  int tspace, keep, t;
  char *s1, *s2, *cpath;
  int fieldSize[128];
  uint8 *mask;
  FILE *mf;
  OneFile *in, *out;

  if (argc != 4)
    { fprintf(stderr, "Usage: alnfilter <in.1aln> <keep.mask> <out.1aln>\n");
      exit(1);
    }
  in = open_Aln_Read(argv[1], 1, &novl, &tspace, &s1, &s2, &cpath);
  if (in == NULL)
    { fprintf(stderr, "alnfilter: cannot open %s\n", argv[1]);
      exit(1);
    }
  mask = (uint8 *) malloc(novl > 0 ? novl : 1);
  mf = fopen(argv[2], "rb");
  if (mf == NULL || mask == NULL || (int64) fread(mask, 1, novl, mf) != novl || fgetc(mf) != EOF)
    { fprintf(stderr, "alnfilter: keep mask %s must hold exactly %lld bytes\n", argv[2], (long long) novl);
      exit(1);
    }
  fclose(mf);
  out = oneFileOpenWriteFrom(argv[3], in, true, 1);
  if (out == NULL)
    { fprintf(stderr, "alnfilter: cannot write %s\n", argv[3]);
      exit(1);
    }
  oneAddProvenance(out, "alnfilter", "0.1", "alnfilter");
  for (t = 0; t < 128; t++)
    fieldSize[t] = in->info[t] != NULL ? in->info[t]->nField * sizeof(OneField) : 0;

  /* header lines up to the first alignment */
  oneGoto(in, 'A', 0);
  while (oneReadLine(in) && in->lineType != 'A')
    { memcpy(out->field, in->field, fieldSize[(int) in->lineType]);
      oneWriteLine(out, in->lineType, oneLen(in), oneString(in));
    }

  /* alignments: the current line is the first 'A' */
  kept = 0;
  for (i = 0; in->lineType == 'A'; i++)
    { keep = (i < novl && mask[i]);
      if (keep)
        { memcpy(out->field, in->field, fieldSize['A']);
          oneWriteLine(out, 'A', oneLen(in), oneString(in));
          kept++;
        }
      while (oneReadLine(in) && in->lineType != 'A')
        if (keep)
          { memcpy(out->field, in->field, fieldSize[(int) in->lineType]);
            oneWriteLine(out, in->lineType, oneLen(in), oneString(in));
          }
      if (in->lineType != 'A')
        break;
    }
  if (i + 1 != novl)
    { fprintf(stderr, "alnfilter: read %lld alignments, header says %lld\n", (long long) (i + 1), (long long) novl);
      exit(1);
    }
  fprintf(stderr, "alnfilter: kept %lld of %lld alignments\n", (long long) kept, (long long) novl);
  oneFileClose(out);
  oneFileClose(in);
  free(mask);
  exit(0);
}
