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
            with open("out.txt", "w") as fo:
                fo.write("%s\n%.6f" % (ans_name, ans_gc * 100.0))
                fo.close()
            f.close()
          
          
          
          
          
          
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
