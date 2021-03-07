#!/usr/bin/env python

#########################################################################################################
## This program will take fasta file of genes (where exons are capitalized and introns are lower case) ##
## and a file of motifs (each motif on it's own line) and will produce a visualiztion of all           ##
## motif locations on each gene using object oriented programming.                                     ##
##########################################################################################################

import cairo
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser("a program to produce two-line fastas")
    parser.add_argument("-f", "--fasta", type=str, help="fasta file of genes to pass through program", required=True)
    parser.add_argument("-m", "--motifs", type=str, help="motif file", required=True)
    parser.add_argument("-o", "--output", type=str, help="name of output (optional), default output uses name of fasta file", required=False)

    return parser.parse_args()
args = get_args()

###assign arguments to variables inside of program
fasta = args.fasta
motifs = args.motifs
out = args.output

#split input name for to get prefix to add to svg output
outputname = re.split('\.', fasta)

#motif disambiguation dictionary
regex_dict = {
    "A":"[Aa]",
    "C":"[Cc]",
    "G":"[Gg]",
    "T":"[TtUu]",
    "U":"[UuTt]",
    "W":"[AaTtUu]",
    "S":"[CcGg]",
    "M":"[AaCc]",
    "K":"[GgTtUu]",
    "R":"[AaGg]",
    "Y":"[CcTtUu]",
    "B":"[CcGgTtUu]",
    "D":"[AaGgTtUu]",
    "H":"[AaCcTtUu]",
    "V":"[AaCcGg]",
    "N":"[AaCcGgTtUu]",
    "Z":"[-]",
}

#############
# Functions #
#############

def find_exon(seq):
    '''this function will return the start and end positions of the exon'''
    exon_tuple = re.search('([A-Z]+)', seq)
    exon = exon_tuple.span()

    return exon

def get_regex(motif):
    '''this function will take a motif and return a regex to find that motif'''
    motif_regex = ""
    for char in motif:
        motif_regex += regex_dict[char]

    return(motif_regex)

def find_positions(seq, reggie):
    '''   '''
    pos_list = []
    motifs = re.finditer(reggie, seq)
    for match in motifs:
        if match is None:
            continue
        else:
            coords = match.span()
            pos_list.append(coords)

    return pos_list

###########
# Classes #
###########

class FastaHeader:
    '''A FastaHeader object stores its start, text, and maybe font.'''
    def __init__(self, context, header, gene_count):
        #data
        self.context = context
        self.header = header
        self.gene_count = gene_count

        #methods
    def draw(self, context):
        context.set_source_rgba(0, 0, 0, 1)
        context.move_to(15, int(self.gene_count)*VERT_PAD + 10)
        context.show_text(self.header)

class Gene:
    '''Gene class – A Gene object should be given enough information that it can figure out how to draw itself.'''
    def __init__(self, context, gene, gene_count):
        #data
        self.context = context
        self.gene = gene
        self.gene_count = gene_count

        #methods
    def draw(self, context):
        '''this method will draw the gene as a horizontal black line, proportional to it's length in nucleotides'''
        context.set_line_width(3)
        context.set_source_rgba(0, 0, 0, 1)
        context.move_to(15, int(self.gene_count)*VERT_PAD + 30)
        context.line_to(15+len(self.gene), int(self.gene_count)*VERT_PAD + 30)
        context.stroke()

class Exon:
    '''An Exon object stores its start, length, and maybe width, similar to Gene.'''
    def __init__(self, context, coords, gene_count):
    #data
        self.coords = coords
        self.gene_count = gene_count

    #methods
    def draw(self, context):
        '''this method will draw the exon on the gene'''
        context.set_source_rgba(0.23,0.25,0.25, 1)
        context.rectangle(15+int(self.coords[0]),int(self.gene_count)*VERT_PAD+20,int(self.coords[1])-int(self.coords[0]),15)        #(x0,y0,x1,y1)
        context.fill()

class Motifs:  
    '''A Motifs object stores a list of motif positions along a given gene and can draw them.'''
    def __init__(self, context, coords_list, R, G, B, gene_count):
    #data
        self.context = context
        self.coords_list = coords_list
        self.gene_count = gene_count
        self.R = R
        self.G = G
        self.B = B
    
    #methods
    def draw(self, context):
        '''this methods will draw rectangles for the motif positions on the gene-line for the coordinates produced by find_positions()'''
        for i in self.coords_list:
            coords = list(i)
            context.set_source_rgba(float(self.R),float(self.G),float(self.B),.7)
            context.rectangle(15+int(coords[0]),int(self.gene_count)*VERT_PAD+20,int(coords[1])-int(coords[0]),15)        #(x0,y0,x1,y1)
            context.fill()

class GeneGroup:
    '''A gene group object manages and organizes the information (gene, exon, and motifs) for each fasta record'''
    def __init__(self):
        
    #data
        self.gene_count = None
        self.header = None
        self.gene = None
        self.exon = None
        self.motifs = []
        self.context = None
    
    #methods
    def draw(self, context):
        self.header.draw(context)
        self.gene.draw(context)
        self.exon.draw(context)
        for i in self.motifs:
            i.draw(context)



#############
# Algorithm #
#############


## parse file to get number of genes and length of longest gene to format context ##
gene_count = 0
gene = ''
seqs = []

#iterate through the file to get the gene count and longest gene for svg dimensions
with open (fasta, "r") as fh:
    for line in fh:
        if line[0] != '>':
            seq = line.strip()
            gene += seq
        else:
            gene_count += 1
            seqs.append(gene)
            gene = ''
        seqs.append(gene)

longest_gene = (len(max(seqs, key=len)))


## figure set up using gene count ##
VERT_PAD = 50 #this number dictates the space between genes, or the 'vertical padding'
WIDTH = int(longest_gene) + 30 #width of figure
HEIGHT = int(gene_count+10) * VERT_PAD #height of figure
surface = cairo.SVGSurface(outputname[0]+'.svg', WIDTH, HEIGHT) #coordinates to display graphic and output name
context = cairo.Context(surface) #create the coordinates you will be drawing on 

#establish motif colors
Rs = ("0.27","1","0.20","0.75","1")
Gs = ("0.5",".63","0.11","0.13",".72")
Bs = ("0.08",".15","1.00","0.06",".83")

## parse through file again to store data as objects and draw the figure ##
motifs_list = []
GENE_COUNT = 0
gene = ''
GENE_GROUPS = []

with open (fasta, "r") as fh, open (motifs, "r") as mt:
    #extract motifs from file into list
    for line in mt:
        motif = line.strip()
        motif = motif.upper()
        if motif not in motifs_list:
            motifs_list.append(motif)
    #start parsing the fasta
    for line in fh:
        if line[0] == '>':
            if gene != '':
                # gene #
                group_obj.gene = Gene(context, gene, GENE_COUNT)
                # exon #
                coords = find_exon(gene)
                group_obj.exon = Exon(context, coords, GENE_COUNT)
                # motifs #
                itR = iter(Rs)
                itG = iter(Gs)
                itB = iter(Bs)
                for motif in motifs_list:
                    #call motif color
                    R = next(itR)
                    G = next(itG)
                    B = next(itB)
                    reggie = get_regex(motif)
                    coords_list = find_positions(gene, reggie)
                    group_obj.motifs.append(Motifs(context, coords_list, R, G, B, GENE_COUNT))
                group_obj.gene_count = GENE_COUNT
                # add gene group to list
                GENE_GROUPS.append(group_obj)
                # reset header and gene seq
                group_obj = GeneGroup()
                header = line.strip()
                GENE_COUNT += 1
                group_obj.header = FastaHeader(context, header, GENE_COUNT)
                gene = ''
            group_obj = GeneGroup()
            header = line.strip()
            #group_obj.gene = Gene(context, gene, GENE_COUNT)
            group_obj.header = FastaHeader(context, header, GENE_COUNT)
            #GENE_GROUPS.append(group_obj)
        else:
            seq = line.strip()
            gene += seq
    #repeat for last sequence
    # gene #
    group_obj.gene = Gene(context, gene, GENE_COUNT)
    coords = find_exon(gene)
    group_obj.exon = Exon(context, coords, GENE_COUNT)
    # motifs #
    itR = iter(Rs)
    itG = iter(Gs)
    itB = iter(Bs)
    for motif in motifs_list:
        #call motif color
        R = next(itR)
        G = next(itG)
        B = next(itB)
        reggie = get_regex(motif)
        coords_list = find_positions(gene, reggie)
        group_obj.motifs.append(Motifs(context, coords_list, R, G, B, GENE_COUNT))
        GENE_GROUPS.append(group_obj)



for gene_group in GENE_GROUPS:
    gene_group.draw(context,)

#draw legend
#make motif colors iterable one last time
itR = iter(Rs)
itG = iter(Gs)
itB = iter(Bs)
motif_counter = 0
#write "legend"
context.set_source_rgba(0, 0, 0, 1)
context.move_to(15,(int(GENE_COUNT)+1)*VERT_PAD+15)
context.show_text("Legend")
for i in motifs_list:
    #call motif color
    R = next(itR)
    G = next(itG)
    B = next(itB)
    #draw legend
    context.set_source_rgba(float(R), float(G), float(B), .9)
    context.rectangle(15,((GENE_COUNT+1)*VERT_PAD)+int(motif_counter+1)*20,40,10)
    context.fill()
    #write motif
    context.move_to(60,((GENE_COUNT+1)*VERT_PAD)+int((motif_counter+1)*20)+10)
    context.show_text(i)
    motif_counter += 1

surface.finish() #close svg file