# Homework 2

EXERCISE 1 DNA
         
          def countNuc(dna):
              counts = {'A': 0, 'C':0, 'G':0,'T':0}
              for nuc in dna:
                  counts[nuc] += 1
              return counts

                      DNAstring = "CCATTCCGTCCGGCCCTCAGCGTGAAACTCCACCTATACGCAAATATCCGCCGAATATAGCGATAAGACCAGTTCGGAGTAATCACTAGAAAAGGACGGCGTGTATCCTCGAGATGCAGATACCAGGTTCATCACTCTAGAAGACCCCATTCGTATGGCAGCCGTATACGCTGTTGCGTGATGCGCGTGCACAGATTTATCCGTACTTGCGGGATTGTCCCGGACCCCACTCCATCCATTCCTCGAAGCCAGGCAGACCTAGCGGTAGCGGGAGGATAGCCCCTAATGGCCCCTCTGTTTGCGCCCGATTAAGGGTCAGATCATTCTAGGTCAGAGAAAGTACATACCCTGTCCCTATACTTTCTACAGCCCCAGCCTGAAAGAATACCAGTAGACTTACTCGTATCCTTCCGAAGTTATGCTCACAGGGCAGGCTTGTGGGTACCACGAGCAGATTGTCGGCATATTATCACCCGCGCACTGATCCACCACGCACGCTCGTGTTGCGATCAATTGTTATAGTGGACCCATGTTGTATCCTAGCTCATCAATAACTTTGGTAAGCGCCAGAAGAGTGCATAGAACGCGGAATCACGATGAGATGCCAATCTACGATCTAGGTGGCACAATCTGGCCCTAACTGAACGACGCTGCTAACTTTATAAATTATCGTAGAGGGATAGACCTTACTAGTGGCTGGAAATACAATGTCTTTTGCGCTCTGTTCCACTGAATCCTAGTCCGGCCCACCACGCAATAATACGGACATAGGCGGAGACAGGGCGGTACCCGCGTTCGTACGTCTACCACGGTTGTAGGAGTCCTCGACTGGGATCCTCTGGGTCGGAGAGAAATGCTCTTCATCTTCCGTTCCATGTACATTGGCTTCATAGTATGCCCGCCGAGTTGCCCAGATCCGTCGTCTACCCTAACAACACAGCGGAGTTCGAATGTAGAAGTCCTTGA"

          result = countNuc(DNAstring)
          print(' '.join([str(val) for key, val in result.items()]))
          
          
          
          
          
          
EXERCISE 2 RNA


          def DNAtoRNA(seq):
                  return seq.replace("T", "U")
          print(DNAtoRNA("CAGTCACGAGAACACGTTCGCAGAACTAGGCGGTCGCTCAGCGCAACGACGCGACTTTCTCACAGGTTCCACTGTCCCCAGTGAGCTCCGAAAAAACCCTCCGCGTTGTCGATACGCCAACTGGTGCCGTCTCCGCCCATATCGGTATATTCTGTCCCGTACATGTAACATTTCTACACTCAACTCATGACTAGATCAGCGACTGACGCAGAGTTTACGTAGTAGACAGACCCAGCGGTAGGCTTTTCCTGTATCATGCATAACGTAGCCGCCATCCAGCGTAAAAATAAGCGGCAATTGCACAGACATCTCTGTTGATTATTTTGGAACGTTGCAATGGAACAGTAGCGCTTAATGTGAGTCACATGACTGTTCGTGACAGCACGGCGTGCCTGTTACTTATAGCAGGGCATCATACAGTCTTAGCGTCGTGTCCATCATTAACACCCTACTGCGTCTATGTCATTAGTAGCAGGATTCGTGCTAAGAGGGTCGTAAAGGCGTTTGTCTTTATAACCTTTTATCTTGTTGAAGGTACAGAATCGGGTAAGTCTATATCGTCACTCTTGGAATTAAAGCAGCATTACTGAGAGAAGCGGACATATTCCGGCACGAGTGCTCTAAACGCTTGGAAGAGTATATGGTCGGGGACACTCGACTATTCAACTGTTCCCCAACCAGACGAGTTGTTCTTTGCAATTCGATGGTCTCATGAAGTTGTCGAGGAAATTATTGTTCTCAACCAAGGGACACACAGCCGGGTCATAAGCCACTAAAAAACATTGACGATGTTCTCCTAAGTGATATCTCATTAAGCGGAAGTATGGGATATTGGCTCAGATATCTAGTGAACTACCTCATATAGCCTCCCCCATTCTGTTCGAACTAGAGTCGTAGCTGTACGAGTCTCACATAGCTATTAGCCGCATCGAGATCAGCTGGGAACCGCCGTCTTTATAACAGGTCTGG")
          
          
 EXERCISE 3 REVC
 
 
          DNA_reversecomplement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

          def reverse_complement(seq):
              return ''.join([DNA_reversecomplement[nuc] for nuc in seq])[::-1]
          print(reverse_complement("ATAGTCCCGCTACGGCGCCTGGGACCCATCGTAGACAGAGATGGAAAAGATGCCGGTCTACAATGAGCGGGTACCGTGGTGGCCGTTATCGCGCCCGTACAAGGACGGTCATTGCAGAGCGCCATACACTCCTGATGTGGCCGCGCCCATGATTGTATAGGTTTATTGCGGACTTGAATGGCGAAAAAGATGTCAACTTGGACGGCGATTGTCCGGCCTACAAGTATCACAGACGCCCCGGGGAGAGCTGTAAGCGAAGGGGATATCTCGCCTTCAATTCCTTCAGATCAAACCGTGTAGATAAAGATTGAAATTCTATGATGGGACATCCATGCCGAGACGGCAAGAGCATGGGTGATTCAGGCGCTCTTGCCCTAGTCGGAAGTTCTCCATCGTGTCAAAAGTTTAGGATCTGGCAGACCGTTCCTCTCCCAATCCCCCCTGGCTTATATTCTAGCGCACGGAACGGGTGTCCTTTATCATTTGCCGTCTACACGCCCACAGCTATCTGACCACCTAGTATCAAAAGTGCCATACTTGCCGTAAATCCTTAGCCCTCCGCTGACTTACATTCCATCTGGGCACGAATTCCCCACTATTCAAACACGAGGCGCGGTTTGCCTACCTCAACGCAGGTGGATTATGCACAATATTCAGTATACCTGAAAGGCCACGCGGAACTTCTGGGCATCACGCGTGTTCCTTACTGATTCTCGTGCCGGTGCCCTCCCCGACAGGGAAATCCTACACGGATTATGGCACGGTGGCCAATCCATCGAGGTTAGTCATATTGTGACAGGGGTTTGGTCGGGGATGTAAGGGAGTCAGAAGTTACGTCGCGTCCACTATACGGAGTTTTTAGCTGTGACCAGATTGAAACTAATGTAAC"))




EXERCISE 4 MENDEL'S FIRST LAW

      import itertools

      def dominant_probability(num_homozygous_dominant, num_heterozygous, num_homozygous_recessive):
          total = num_homozygous_dominant + num_heterozygous + num_homozygous_recessive

          recessive_probability = (num_homozygous_recessive / total) * ((num_homozygous_recessive - 1) / (total - 1))
          heterozygous_probability = (num_heterozygous / total) * ((num_heterozygous - 1) / (total - 1))
          hetero_recessive_probability = (num_heterozygous / total) * (num_homozygous_recessive / (total - 1)) + (num_homozygous_recessive / total) * (num_heterozygous / (total - 1))

          recessive_total = recessive_probability + heterozygous_probability * (.25) + hetero_recessive_probability * (.5)

          return (1 - recessive_total)

      if __name__ == "__main__":

          result = round(dominant_probability(17.0, 19.0, 26.0), 5)

          print (result)
 
          f = open('workfile.txt', 'w')
          f.write(str(result))
          
          
          
Gc content exercise 5
          
        with open("rosalind_gc.txt", "r") as f:
            seq_name = f.readline().replace("\n", "")
            seq_content = ""
            ans_gc = 0.0
            ans_name = ""
            while True:
                line = f.readline().replace("\n", "")
                if (line == "" or line[0] == '>'):
                    now_gc = (seq_content.count('G') + seq_content.count('C'))* 1.0 / len(seq_content)
                    if (now_gc > ans_gc):
                       ans_gc = now_gc
                       ans_name = seq_name
                    if (line == ""):
                       break
                    seq_name = line[1:]
                    seq_content = ""
                    continue
                seq_content = seq_content + line
            print(ans_name, ans_gc * 100.0)
            
          
          
          
          
          
          
EXERCISE 6 SUBS
          
          
      s="CGTAACCCCATGGCAAAACCCCAGCAAACCCCAAGAGCAACCCCACAAACCCCAATAACCCCAAACCCCAAACCCCAGGCATTCACGTAACCCCACCATAACCCCAAACCCCAGAACCCCAAACCCCATGAACCCCACGTTCAACCCCAAACCCCACAGGCAACCCCAAAACCCCAAACCCCACAACCCCAGTAACCCCAAGAACCCCATAACCCCAGGGACAACCCCATTGAACCCCAGGAACCCCACAACCCCAGAAACCCCACGGAGATTACAACCCCACGCTCTAAACCCCATGTCAAACCCCAAACCCCAAACCCCAAACCCCACAAACCCCAAACCCCAAACCCCAGGCAACCCCAGGGCTACTGCTCAACCCCATGAAACCCCACAAACCCCATTAGCTGAAACCCCAGAACCCCATAACCCCAGTCAGCTGAACCCCAGGCAGGTTTCAACCCCAAACCCCATCAAACCCCAAAACCCCATAACCCCACATAACCCCAATTCGAACCCCAAAACCCCAAACCCCAAAACCCCAGACCTAACCCCACGCTCTGAACCCCATGATGTAACCCCAAACCCCAGAACCCCACTGCGCTACACAACCCCAGCAACCCCAACCAACCCCAAACCCCATGTTCACGAACCCCAGAACCCCAGGCCAACCCCACAAACCCCAGGTCAGGGACCAAACCCCAAACCCCATTTCATGAACCCCAAGCTAACCCCAGGAACCCCAGTCAAACCCCATTTAACCCCACAACCCCAGCGAAGAGAACCCCACAACCCCAAACCCCAAACCCCAAACCCCAGGAACCCCACGGGGAACCCCAAACCCCAGGCAACCCCAAACCCCAATAACCCCAAACCCCATATTAACCCCA"
      t = "AACCCCAAA"

      for position in range(len(s)):
          if s[position:].startswith(t):
              print(position+1)
              
              
              
              
              
 Exercise 7 LIA
              
              
      from math import factorial as FACT
      k = 5
      N = 8
      P = 2**k
      probability = 0
      for i in range(N,P+1):
          prob = (FACT(P)/(FACT(i)* FACT(P-i)))*(0.25**i)*(0.75**(P-i))
          probability += prob
      print(round(probability,3))
      
      
  EXERCISE 8 IEV 
      
      
from enum import Enum
class Allele(Enum):
    DOMINANT = 1,
    RECESSIVE = 2,
    HET = 3


def mate(x,y):
    if x == Allele.DOMINANT and y == Allele.DOMINANT:
        return [Allele.DOMINANT]*4
    if x == Allele.RECESSIVE and y == Allele.RECESSIVE:
        return [Allele.RECESSIVE]*4
    if x == Allele.HET and y == Allele.HET:
        return [Allele.DOMINANT, Allele.HET, Allele.HET, Allele.RECESSIVE]
    if x == Allele.HET:
        return [Allele.HET, Allele.HET, y, y]
    if y == Allele.HET:
        return [Allele.HET, Allele.HET, x, x]
    else:#DOMINANT AND RECISSIVE
        return [Allele.HET]*4


def get_pairs(organisms):
    if len(organisms) < 2:
        return list()
    first = organisms[0]
    pairs = [(first, organisms[second_index]) for second_index in range(1,len(organisms))]
    other_pairs = get_pairs(organisms[1:])
    pairs.extend(other_pairs)
    return pairs

def f(AAAA, AAAa, AAaa, AaAa, Aaaa, aaaa):
    children = list()
    for i in range(AAAA): 
        children.extend(mate(Allele.DOMINANT, Allele.DOMINANT))
    for i in range(AAAa):
        children.extend(mate(Allele.DOMINANT, Allele.HET))
    for i in range(AAaa):
        children.extend(mate(Allele.DOMINANT, Allele.RECESSIVE))
    for i in range(AaAa):
        children.extend(mate(Allele.HET, Allele.HET))
    for i in range(Aaaa):
        children.extend(mate(Allele.HET, Allele.RECESSIVE))
    for i in range(aaaa):
        children.extend(mate(Allele.RECESSIVE, Allele.RECESSIVE))

    count_dominant_presence = 0
    for child in children:
        if child is not Allele.RECESSIVE:
            count_dominant_presence += 1

    return count_dominant_presence/2
AAAA = 19047  
AAAa = 16432
AAaa = 17007
AaAa = 18621
Aaaa = 16638
aaaa = 16427

print(f(AAAA, AAAa, AAaa, AaAa, Aaaa, aaaa))
      
      
exercise 9 consensus and profile
      
      
      def read_fasta(fp):
          name, seq = None, []
          for line in fp:
              line = line.rstrip()
              if line.startswith(">"):
                  if name: yield (name, ''.join(seq))
                  name, seq = line, []
              else:
                  seq.append(line)
          if name: yield (name, ''.join(seq))

      data_list = []
      with open('rosalind_cons.txt') as fp: 
          for name, seq in read_fasta(fp):
              data_list.append(seq)
      length = len(data_list)
 
      L = len(seq)
 
      P = [[0 for x in range(L)] for y in range(4)] 
 
      Q = ['A','C','G','T']
 
      for x in range(L):
          for y in range(4):
              for z in range(length):
                  P[y][x] = P[y][x] + data_list[z][x].count(Q[y])
 
      domi = ""
      for x in range(L):
          MAX = 0
          for y in range(4):  
              if P[y][x]>=P[MAX][x]:
                  MAX = y       
             
          if MAX == 0:
              domi = domi+'A'
          elif MAX ==1:
              domi = domi+'C'
          elif MAX ==2:
              domi = domi+'G'
          elif MAX ==3:
              domi = domi+'T'

      print('%s\n%s\n%s\n%s\n%s' %(
          domi,'A: '+str(P[0]).strip('[]').replace(',',''),'C: '
          +str(P[1]).strip('[]').replace(',',''),'G: '
          +str(P[2]).strip('[]').replace(',',''),'T: '
          +str(P[3]).strip('[]').replace(',','')))




Exercise 10 prob



     import math

     def RandomString(strRandomString, stringArray):
         strRandomString = strRandomString.upper()
         cg = len(strRandomString.replace('A', '').replace('T', ''))
         at = len(strRandomString.replace('C', '').replace('G', ''))
         inputArray = stringArray.split()
         outputArray = []
         for i in range(0, len(inputArray)):
             prob = cg * math.log10(float(inputArray[i]) / 2) + at * math.log10((1 - float(inputArray[i])) / 2)
             outputArray.append(round(prob, 3))
         return outputArray




     print(' '.join(map(str,  RandomString('CATGGTGTACACGACCTCTCGTTTTCCTTACCAAACCGATCAGAGCCCTCATGTACGCGCTGGAGACGAATAAGCAGAGTTACTGTAATTACTACTAAA', '0.062 0.130 0.190 0.238 0.266 0.341 0.371 0.439 0.504 0.573 0.590 0.680 0.733 0.759 0.804 0.854 0.920'))))
