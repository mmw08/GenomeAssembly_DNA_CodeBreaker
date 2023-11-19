#resources:
#https://github.com/navjindervirdee/Genome-Sequencing/blob/master/Genome%20Assembly%20using%20De%20Bruijn%20Graphs%20from%20K-Mer%20Composition/GenomeAssembly.java
#https://github.com/xuwd11/Coursera-Bioinformatics/blob/master/61_02_TrieMatching.py

import sys
from collections import defaultdict
import time

reads_file= sys.argv[1]
reference_file= sys.argv[2]
Mismatch_Flag=sys.argv[3] #Should only be Y/N


def DNA_CODE(reads_file, reference_file, Mismatch_Flag,Name_of_output_file):
    #If the reads is exact reads without missmatch
    if Mismatch_Flag.upper()=="N":
        
        #Function to of  construct the the debruijnGraph where the length of kemer is length of first element of reads
        def deBruijnGraph(kmers):
            graph = {}
            for i in range(len(kmers)):
                try:
                    graph[kmers[i][:-1]].append(kmers[i][1:])
                except:
                    graph[kmers[i][:-1]] = [kmers[i][1:]]
            return graph

        # Read the kmers  from the 'reads.txt' file
        with open(reads_file, 'r') as file:
            kmers = [line.strip() for line in file]

         # build the graph where node is the length of read
        graph = deBruijnGraph(kmers) 

        input_edges = []
        # find the edge btw the nodes and add them to the list of edges
        for key, values in graph.items():
            for v in values:
                input_edges.append(f"{key} -> {v}")

        #counting the number of outgoing edges (out-degree) and the number of incoming edges (in-degree). The resulting degrees dictionary will have the node identifiers as keys, and the values will represent the difference between out-degrees and in-degrees for each node.
        degrees = defaultdict(int)
        for k in graph:
            for v in graph[k]:
                degrees[k] += 1
                degrees[v] -= 1
        #identifying the source and sink nodes in a directed graph based on the calculated in-degrees and out-degrees stored in the degrees dictionary
        source = [k for k, v in degrees.items() if v == 1][0]
        sinc = [k for k, v in degrees.items() if v == -1][0]

        list(graph)
        start = list(graph)[0]

        #ensures that the sink node is connected to the source node in the graph, so that the traversal can form a cycle
        if sinc in graph.keys():
            graph[sinc].append(source)
        else:
            graph[sinc] = [source]

        #implementing a traversal to find cycles in a  graph, particularly for Eulerian cycles
        cycles = {}
        while graph:
            current = next(iter(graph))
            cycle = [current]
            cycles[current] = cycle
            while current in graph:
                _next = graph[current][0]
                del graph[current][0]
                if len(graph[current]) == 0:
                    del graph[current]
                current = _next
                cycle.append(_next)
        #performs a traversal of the tree represented by the cycles dictionary starting from the specified root node to form a cycle list contains the nodes visited during the traversal
        def traverse(tree, root):
            out = []
            for r in tree[root]:
                if r != root and r in tree:
                    out += traverse(tree, r)
                else:
                    out.append(r)
            return out

        cycle = traverse(cycles, start)
        for i in range(1, len(cycle)):
            if cycle[i-1] == sinc and cycle[i] == source:
                boarder = i
        path = cycle[boarder:]+cycle[1:boarder]

        # Assemble the genome sequence
        assembled_genome = ''.join([s[0] for s in path]) + ''.join(path[-1][1:])

        # Write the assembled genome to the 'print.txt' file
        with open(Name_of_output_file, 'w') as output_file:
            output_file.write(assembled_genome)

            #check if the user has misreads
    elif Mismatch_Flag.upper()=="Y" :
        #read the reads file 
        with open(reads_file, 'r') as f:

            lines= f.readlines()
            reads=[]
            for i in lines:
                reads.append(i[:-1])
        Reference_genome=""

        #read the reference file
        with open(reference_file, 'r') as f:
            lines= f.readlines()
            
            for i in lines:
                Reference_genome+=i[:-1]

        #FUNCTION THAT CONSTRUCT THE TRIE usig the reads file and return dictiony that store the nodes
        def trieconstructor(read):
            trie = {}
            trie[0] = {}
            node_position = 1
            for read_i in read:
                current = 0
                for j in range(len(read_i)):
                    p = read_i[j]
                    if p in trie[current]:
                        current = trie[current][p]
                    else:
                        trie[node_position] = {}
                        trie[current][p] = node_position
                        current = node_position
                        node_position += 1
                trie[current]={'$':read_i}
            return trie

        #check if the part of the ref exist in the read 
       

        def check_sub_m(s, trie, allow_mismatch=True):
            node = 0
            i = 0
            sub = ""
            mismatch_occurred = not allow_mismatch
            
            while True:
                if '$' in trie[node]:  # Found a pattern in the trie
                    return sub
                
                if i < len(s):
                    symbol = s[i]
                    if symbol in trie[node]:  # Exact match
                        sub += symbol
                        node = trie[node][symbol]
                        i += 1
                    elif not mismatch_occurred:  # Allow for one mismatch
                        mismatch_occurred = True
                        sub += symbol
                        i += 1
                    else:  # Mismatch has already occurred, no more allowed
                        return None
                else:
                    return None
        #match reads and string 
        def triematching(s,pattern):
            
            trie = trieconstructor(pattern)
            
            dict_postion = {}
            n = len(s)
            for i in range(n):
                check_i = check_sub_m(s[i:], trie)
                if check_i is not None:
                    if check_i in dict_postion:
                        dict_postion[check_i].append(i)
                    else:
                        dict_postion[check_i]=[i]
            return dict_postion

        result = triematching(Reference_genome, reads)
        
        #build a string of length of reference  and return a list contains the postion for each pattern 
        Final_re=  '-' * len(Reference_genome)

        Final_re_list = list(Final_re)

        for replacement, positions in result.items():
            for pos in positions:
                Final_re_list[pos:pos+len(replacement)] = replacement

        # Convert the list back to a string
        Final_re = ''.join(Final_re_list)    
        
        #write the string in the the file and if there is  multiple contig we wirte the in different lines
        with open(Name_of_output_file,"w") as f:
            for i in range(len(Final_re)):
                if i!=0 and Final_re[i]=="-" and Final_re[i-1]in ["A","C","T","G"]:
                    f.write("\n")
                elif Final_re[i]=="-":
                    continue
                else:

                    f.write(Final_re[i])
       

    else:
        print("Misreads only accept N/Y")    


DNA_CODE(reads_file,reference_file,Mismatch_Flag,"assembledreads.txt")
