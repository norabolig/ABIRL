import os
import argparse

def main():
    parser= argparse.ArgumentParser(description="Loop through files to align based on a reference")
    parser.add_argument("reference",help="Reference fits file")
    parser.add_argument("list",help="list file")
    parser.add_argument("prefix",help="output fits file prefix")
    parser.add_argument("--istart",type=int,default=0,help="starting vaue for iterator and output file")
    parser.add_argument("--catalogue", help="Optional: CSV filename to save star matches")
    parser.add_argument("--thresh", type=float, default=0.05, help="Detection threshold (lower = more sensitive)")
    parser.add_argument("--max_src", type=int, default=200, help="Max stars to detect")
    parser.add_argument("--fill", type=float, default=0.0, help="Value for empty pixels")
    parser.add_argument("--wcs", type=bool, default=False, help="Align by WCS")

    args = parser.parse_args()

    listfile=args.list

    fh=open(listfile,"r")

    iter = args.istart*1
    for line in fh:
        shiftfile=line.rstrip()

        if args.wcs: os.system("python3 ~/ABIRL/imalign_adv.py {} {} {} --wcs --catalogue {} --thresh {} --max_src {} --fill {}".format(args.reference,shiftfile,args.prefix+str(iter).zfill(6)+".fits",args.catalogue,args.thresh,args.max_src,args.fill))
        else: os.system("python3 ~/ABIRL/imalign_adv.py {} {} {} --catalogue {} --thresh {} --max_src {} --fill {}".format(args.reference,shiftfile,args.prefix+str(iter).zfill(6)+".fits",args.catalogue,args.thresh,args.max_src,args.fill))
        iter+=1

    fh.close()

if __name__ == "__main__":
    main()
