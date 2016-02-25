import pandas,getopt,sys

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
      print 'test.py -i <inputfile> -o <outputfile>'
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
            print 'test.py -i <inputfile> -o <outputfile>'
            sys.exit()
      elif opt in ("-i", "--ifile"):
            inputfile = arg
      elif opt in ("-o", "--ofile"):
            outputfile = arg
            
    popstats = pandas.read_hdf(inputfile,'popstats')
    print popstats.columns.values
    popstats['evg']=4*popstats['mu']
    g1 = popstats.groupby(['mu','sige','r','sigmu']).mean()
    g1.to_csv(outputfile,sep="\t",header=True)

if __name__ == "__main__":
   main(sys.argv[1:])
