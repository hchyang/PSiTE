#!/usr/bin/env python

#########################################################################
# Author: Hechuan
# Created Time: 2017-04-04 18:00:34
# File Name: simulate_somatic_vars.py
# Description: 
#########################################################################

import os
import sys
import re
import pickle
import numpy
import argparse
import copy
import logging
import trunk_vars

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

class Tree:
    count=0
    cnv_c2=0
    def __init__(self,name=None,lens=None,left=None,right=None,top=None,snvs=None,accumulated_snvs=None,cnvs=None,accumulated_dels=None,C='0.0.0'):
        self.name=name
        self.lens=lens
        self.left=left
        self.right=right
        self.top=top
        self.snvs=snvs                         #it contains position of each snv that occured on its top branch
        self.accumulated_snvs=accumulated_snvs #it contains pos for all snvs on the lineage leading to that node 
#it's a list of dictionary and each dictionary contains those keys {start,end,copy,leaves_count,pre_snvs,new_copies} of each cnv that occured on its top branch
        self.cnvs=cnvs                         
        self.accumulated_dels=accumulated_dels #it contains [start,end] for each del on the lineage leading to that node 
        self.C=C
        Tree.count+=1

    def add_node(self,node):
        '''
        This method is used when constructing a tree.
        After adding node to the growing tree, return the node added as the new tree.
        '''
        if self.left==None:
            self.left=node
            node.top=self
            self=node
        elif self.right==None:
            self.right=node
            node.top=self
            self=node
        else:
            print('Should not be here!1')
        return self

####################################################################################################

    def add_snv_cnv(self,start=0,end=1,inherent_snvs=[],inherent_dels=[],inherent_cnvs=[],
                    snv_rate=None,cnv_rate=None,del_prob=None,cnv_length_lambda=None,cnv_length_max=None,copy_max=None):
        '''
        Randomly put SNVs and CNVs on a phylogenetic tree.
        For amplifications, we will build a new tree for every new copy and use this method to add SNVs/CNVs on the new tree.
        NOTE: 1. snv_rate/cnv_rate are correlated with start end. Make sure start-end==1.
              2. We assume for each CNVs, its start is inclusive, but its end is not.
        '''
#rescale mutation rate according the length of the sequence
        length=end-start
        logging.debug('New node with length: %s',self.lens)
        logging.debug('Structure: %s',self.tree2newick())
        self.snvs=[]
        self.cnvs=[]
        self.accumulated_snvs=[]
        self.accumulated_dels=[]

        if self.top == None: 
#root node, may have inherent_snvs
            self.snvs=inherent_snvs[:]
            self.cnvs=inherent_cnvs[:]
            self.accumulated_snvs=inherent_snvs[:]
            self.accumulated_dels=inherent_dels[:]
        else:
#non-root node inherits snvs/cnvs from its top nodes 
            if self.top.accumulated_snvs != None:
                self.accumulated_snvs=self.top.accumulated_snvs[:] 
            if self.top.accumulated_dels != None:
                self.accumulated_dels=self.top.accumulated_dels[:]

        mutation_rate=(snv_rate+cnv_rate)*length
        if self.lens != None and mutation_rate > 0:
            snv_prob=snv_rate/(snv_rate+cnv_rate)
            mutation_waiting_times=waiting_times(span=self.lens,rate=mutation_rate)
            for waiting_t in mutation_waiting_times:
                pos=numpy.random.uniform(start,end)
                if numpy.random.uniform()<snv_prob:
#snv
                    if self.accumulated_dels:
                        for del_start,del_end in self.accumulated_dels:
                            if not del_start<=pos<del_end:
                                self.snvs.append(pos)
                                self.accumulated_snvs.append(pos)
                    else:
                        self.snvs.append(pos)
                        self.accumulated_snvs.append(pos)
                    logging.debug('New SNV: %s',pos)
                    logging.debug('The length of the tree with new SNV: %s',self.lens)
                    logging.debug('Structure: %s',self.tree2newick())
                else:
#cnvs
#if the new cnv overlap with accumulated_dels, compare it with the accumulated dels 
#and only keep those new regions.
                    cnv_start=pos
                    cnv_length=numpy.random.exponential(cnv_length_lambda)
                    if cnv_length>cnv_length_max:
                        cnv_length=cnv_length_max
                    cnv_end=cnv_start+cnv_length
                    if cnv_end>end:
                        cnv_end=end
                    leaves_count=self.leaves_number()
                    new_cnvs=[[cnv_start,cnv_end]]
                    logging.debug('New CNV: %s',str(new_cnvs))
                    logging.debug('Pre deletions: %s',str(self.accumulated_dels))
                    for cnv in new_cnvs: #We need to modify new_cnvs in place
                        for del_start,del_end in self.accumulated_dels:
                            if cnv[0]<del_start:
                                if del_start<=cnv[1]<=del_end:
                                    cnv[1]=del_start
                                elif cnv[1]>del_end:
                                    new_cnvs.append([del_end,cnv[1]])
                                    cnv[1]=del_start
                            elif del_start<=cnv[0]<=del_end:
                                if del_start<=cnv[1]<=del_end:
                                    new_cnvs.remove(cnv)
                                    break
                                elif cnv[1]>del_end:
                                    cnv[0]=del_end
                                else:
                                    print('Should not be here!3')
                    logging.debug('New CNVs after comparing with pre deletions: %s',str(new_cnvs))
                    if len(new_cnvs)==0 or len(new_cnvs[0])==0:
                        continue
########################################################################################################################
                    if numpy.random.uniform()<del_prob:
#the new cnv is a deletion
                        logging.debug('New CNVs are deletions.')
                        for del_start,del_end in new_cnvs:
#output pre_snvs to self.cnvs, so it can be used to correct the count of snvs 
                            pre_snvs=[]
                            for snv in self.accumulated_snvs:
                                if del_start<=snv<del_end:
                                    self.accumulated_snvs.remove(snv)
                                    if snv in self.snvs:
                                        self.snvs.remove(snv)
                                    else:
                                        pre_snvs.append(snv)
                            cnv={'seg':[start,end],
                                 'start':del_start,
                                 'end':del_end,
                                 'copy':-1,
                                 'leaves_count':leaves_count,
                                 'pre_snvs':pre_snvs,
                                 'new_copies':[]}
                            self.cnvs.append(cnv)
                            self.accumulated_dels.append([del_start,del_end])
                    else:
#the new cnv is an amplification
                        logging.debug('New CNVs are amplifications.')
                        cnv_copy=numpy.random.random_integers(copy_max)
                        for amp_start,amp_end in new_cnvs:
                            amp_length=amp_end-amp_start
#collect the old snvs on cnvs. Those snvs are the snvs on the ancestor lineage leading to segment, and locate in the segment.
                            pre_snvs=[]
                            for snv in self.accumulated_snvs:
                                if amp_start<=snv<amp_end:
                                    pre_snvs.append(snv)

#collect the new snvs on cnvs
                            new_copies=[]
                            for i in range(cnv_copy):
                                segment=Tree(name=self.name,lens=float(self.lens)-waiting_t)
                                if self.left != None:
                                    segment.left=copy.deepcopy(self.left)
                                    segment.left.top=segment
                                if self.right != None:
                                    segment.right=copy.deepcopy(self.right)
                                    segment.right.top=segment
                                new_copies.append(segment)
                            cnv={'seg':[start,end],
                                 'start':amp_start,
                                 'end':amp_end,
                                 'copy':cnv_copy,
                                 'leaves_count':leaves_count,
                                 'pre_snvs':pre_snvs,
                                 'new_copies':new_copies}
                            self.cnvs.append(cnv)
        for cnv in self.cnvs:
            if cnv['copy']>0: 
                for segment in cnv['new_copies']:
                    segment.add_snv_cnv(start=cnv['start'],end=cnv['end'],inherent_snvs=cnv['pre_snvs'],
                                        snv_rate=snv_rate,cnv_rate=cnv_rate,del_prob=del_prob,
                                        cnv_length_lambda=cnv_length_lambda,cnv_length_max=cnv_length_max,copy_max=copy_max)

        if self.left != None:
            self.left.add_snv_cnv(start=start,end=end,inherent_snvs=[],
                                  snv_rate=snv_rate,cnv_rate=cnv_rate,del_prob=del_prob,
                                  cnv_length_lambda=cnv_length_lambda,cnv_length_max=cnv_length_max,copy_max=copy_max)
        if self.right != None:
            self.right.add_snv_cnv(start=start,end=end,inherent_snvs=[],
                                   snv_rate=snv_rate,cnv_rate=cnv_rate,del_prob=del_prob,
                                   cnv_length_lambda=cnv_length_lambda,cnv_length_max=cnv_length_max,copy_max=copy_max)

    def all_cnvs_collect(self):
        '''
        Return a list of all cnvs on the tree.
        '''
        all_cnvs=[]
        if self.cnvs != None:
            all_cnvs=self.cnvs[:]
            tmp_cnvs=self.cnvs[:]
            for cnv in tmp_cnvs: #We should use tmp_cnvs instead of all_cnvs here, as all_cnvs is growing in this 'for' loop.
                if cnv['copy']>0:
                    for cp in cnv['new_copies']:
                        all_cnvs.extend(cp.all_cnvs_collect())
        if self.left != None:
            all_cnvs.extend(self.left.all_cnvs_collect())
        if self.right != None:
            all_cnvs.extend(self.right.all_cnvs_collect())
        return all_cnvs

########################################################################################################################
#self.cnvs has those keys: [start,end,copy,leaves_count,pre_snvs,new_copies] of each cnv that occured on its top branch
#how to count snvs in deletions
    def all_snvs_allele_count(self):
        '''
        It will return a diction of all SNVs in the tree. 
        {pos:alt_allele_count,...}
        We will summary the allele count of all snvs on both main tree and all
        subtrees (new copy of cnv).
        There are three kind of snvs: 
        1) on the main tree 
        2) on the subtree (pre_snvs and new snvs)
        3) in deletions (pre_snvs)
        '''
        all_alt_count={}
        if self.snvs!=None:
            for snv in self.snvs:
                all_alt_count[snv]=self.leaves_number()
        if self.cnvs!=None:
            for cnv in self.cnvs:
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        all_alt_count=merge_two_dict(all_alt_count,cp.all_snvs_allele_count())
                else:  #deletion
                    pre_snvs_dict={}
                    for snv in cnv['pre_snvs']:
                        pre_snvs_dict[snv]=-self.leaves_number()
                    all_alt_count=merge_two_dict(all_alt_count,pre_snvs_dict)

        if self.left!=None:
            all_alt_count=merge_two_dict(all_alt_count,self.left.all_snvs_allele_count())
        if self.right!=None:
            all_alt_count=merge_two_dict(all_alt_count,self.right.all_snvs_allele_count())
        return all_alt_count

    def snvs_alt_count(self):
        '''
        In this function, I will count the leaves of each snv and return the list of sorted snvs:
        [[pos,alt_count],...]
        '''
        snvs_alt_count_dict=self.all_snvs_allele_count()
        snvs=list(snvs_alt_count_dict.keys())
        snvs.sort(key=float)
        snvs_alt_count=[[snv,snvs_alt_count_dict[snv]] for snv in snvs]
        return snvs_alt_count

    def leaves_number(self):
        '''
        After this method, all nodes will have the attribute of leaves count.
        '''
        if not hasattr(self,'leaves_count') or self.leaves_count == None:
            if self.left==None and self.right==None:
                self.leaves_count=1
            else:
                self.leaves_count=0
                if self.left!=None:
                    self.leaves_count+=self.left.leaves_number()
                if self.right!=None:
                    self.leaves_count+=self.right.leaves_number()
        return self.leaves_count

#FIXME: Should we store all the leaves' names in each node?
    def leaves_names(self):
        '''
        After this method, all nodes will have the attribute of leaves names.
        '''
        if not hasattr(self,'leaves_names') or self.leaves_names == None:
            if self.left==None and self.right==None:
                self.leaves_names=[self.name]
            else:
                self.leaves_names=[]
                if self.left!=None:
                    self.leaves_names+=self.left.leaves_names()
                if self.right!=None:
                    self.leaves_names+=self.right.leaves_names()
        return self.names
    
    def snvs_freq_cnvs_profile(self,ploid=None,snv_rate=None,cnv_rate=None,del_prob=None,
                               cnv_length_lambda=None,cnv_length_max=None,copy_max=None,
                               trunk_snvs=None,trunk_dels=None,trunk_cnvs=None):
        '''
        Produce the true frequency of SNVs in the samples.
        It's a warpper for generating SNVs/CNVs on a tree and summarize their frequency.
        '''
        all_cnvs=[]
        all_snvs_alt_counts=[]
        all_snvs_alt_freq=[]
        all_leaves=self.leaves_number()
        background=all_leaves*ploid
        logging.debug('Your tree is: %s',self.tree2newick())
#collect all snvs and cnvs
        for i in range(ploid):
            logging.info(' Simulate tree %s (total: %s)',i+1,ploid)
            hap_tree=copy.deepcopy(self)
            hap_trunk_snvs=[]
            hap_trunk_dels=[]
            hap_trunk_cnvs=[]
            if i in trunk_snvs:
                hap_trunk_snvs=trunk_snvs[i]
            if i in trunk_dels:
                hap_trunk_dels=trunk_dels[i]
            if i in trunk_cnvs:
                hap_trunk_cnvs=trunk_cnvs[i]
            hap_tree.add_snv_cnv(inherent_snvs=hap_trunk_snvs,inherent_dels=hap_trunk_dels,inherent_cnvs=hap_trunk_cnvs,
                                 snv_rate=snv_rate,cnv_rate=cnv_rate,del_prob=del_prob,
                                 cnv_length_lambda=cnv_length_lambda,cnv_length_max=cnv_length_max,copy_max=copy_max)
            all_snvs_alt_counts.extend(hap_tree.snvs_alt_count())
            all_cnvs.extend(hap_tree.all_cnvs_collect())

        all_cnvs.sort(key=lambda cnv: cnv['start'])

#construct depth profile list, assuming the whole region start with 0 and end with 1.
        all_pos_changes=cnvs2break_points(all_cnvs)
        all_pos_changes.insert(0,[0,background])
        all_pos_changes.append([1,-background])

        depth_profile=pos_changes2region_profile(all_pos_changes)
        all_snvs_alt_counts.sort(key=lambda snv: snv[0])
        region_mean_ploid=0
        for snv in all_snvs_alt_counts:
            while snv[0]>=all_pos_changes[0][0]:
                change=all_pos_changes.pop(0)
                region_mean_ploid+=change[1]
            all_snvs_alt_freq.append([snv[0],snv[1]/region_mean_ploid])
        return all_snvs_alt_freq,all_cnvs,depth_profile

####################################################################################################

    def tree2newick(self,lens=False,attr=None):
        '''
        Convert tree structure to string in Newick/NHX format.
        '''
        newick_str=''

        if self.left != None:
            newick_str+='(' + self.left.tree2newick(lens=lens,attr=attr)
        if self.name==None:
            newick_str+=','
        else:
            newick_str+=self.name
        if self.right != None:
            newick_str+=self.right.tree2newick(lens=lens,attr=attr) + ')'

        if self.lens!=None and lens:
            newick_str+= ':' + self.lens
        if attr!=None and lens==True:
            newick_str+='[&&NHX:W=1.0'
            for attribute in attr:
                if getattr(self, attribute):
                    newick_str+=':{}={}'.format(attribute,getattr(self, attribute))
            newick_str+=']'
        return newick_str

    def snvs_with_carriers(self):
        snvs_carriers=[]
        if self.snvs != None:
            for snv in self.snvs:
                snvs_carriers.append([snv]+self.leaves_names())
        return snvs_carriers

    def snvs_all_raw(self):
        snvs_all=self.snvs_with_carriers()
        if self.left != None:
            snvs_all.extend(self.left.snvs_all_raw())
        if self.right != None:
            snvs_all.extend(self.right.snvs_all_raw())
        return snvs_all

################################################################################
#FIXME this method is not currect because of the effect of CNVs
#let's sort snvs by their positions and output genotypes for all samples
#in each locus
    def snvs_all(self):
        snvs_all=self.snvs_all_raw()
        snvs_all.sort(key=lambda snv: snv[0])
        all_leaves=self.leaves_names()
        all_genotypes=[]
#FIXME: try to make sure there will not be two snvs at the same postion 
        for snv in snvs_all:
            genotype=[snv[0]]
            for leaf in all_leaves:
                if leaf in snv[1:]:
                    genotype+=[1]
                else:
                    genotype+=[0]
            all_genotypes+=[genotype]
        return all_genotypes
################################################################################

    def highlight_snvs(self,snvs):
        if self.left != None:
            self.left.highlight_snvs(snvs)
        if self.right != None:
            self.right.highlight_snvs(snvs)
        if self.snvs !=None and set(self.snvs).intersection(snvs):
            self.C='255.0.0'

def waiting_times(span=None,rate=None):
    elapse=0.0
    waiting_times=[]
    span=float(span)
    if rate>0:
        while elapse<span:
            elapse+=numpy.random.exponential(1/rate)
            if elapse<span:
                waiting_times.append(elapse)
    return waiting_times

def merge_two_dict(dict1={},dict2={}):
    new_dict={}
    for key in set.union(set(dict1),set(dict2)):
        if key in dict1 and key in dict2:
            new_dict[key]=dict1[key]+dict2[key]
        elif key in dict1:
            new_dict[key]=dict1[key]
        else:
            new_dict[key]=dict2[key]
    return new_dict
    
def cnvs2break_points(cnvs):
    '''
    Return a list of lists. Each sublist contain two elements. The first is the postion and the second 
    is the copy number CHANGES across all the samples between that positon and the next position.
    [[pos,relative_copy_number_change],...]
    '''
    pos_changes=[]
    for cnv in cnvs:
        change=cnv['copy']*cnv['leaves_count']
        pos_changes.extend([[cnv['start'],change],[cnv['end'],-change]])
    pos_changes.sort(key=lambda pos: pos[0])
    return pos_changes

def pos_changes2region_profile(pos_changes):
    '''
    Convert [[pos,change]...] to [[start,end,current]...]
    '''
    pos_changes_copy=pos_changes[:]
    profile=[]
    current=0
    while pos_changes_copy:
        pos,change=pos_changes_copy.pop(0)
        start=pos
        current+=change
        if pos_changes_copy:
            end=pos_changes_copy[0][0]
            if start!=end:
                profile.append([start,end,current])
    return profile

def newick2tree(newick=None):
    leaf_name_re=re.compile('\w+:')
    lens_re=re.compile(':[0-9.]+')
    while newick != ';':
        if newick.startswith('('):
            if 'mytree' in vars():
                mytree=mytree.add_node(Tree())
            else:
                mytree=Tree()
            newick=newick[1:]
        elif leaf_name_re.match(newick):
            m=leaf_name_re.match(newick)
            index=m.span()
            leaf_name=newick[:index[1]-1]
            newick=newick[index[1]-1:]
            mytree=mytree.add_node(Tree(name=leaf_name))
        elif lens_re.match(newick):
            m=lens_re.match(newick)
            index=m.span()
            lens=newick[1:index[1]]
            newick=newick[index[1]:]
            mytree.lens=lens
        elif newick.startswith((',',')')):
            newick=newick[1:]
            mytree=mytree.top
        else:
            print('Should not be here!4')
            print('Check your newick tree! Or maybe there are something I do not know about newick!')
    return mytree
    
def simulate_sequence_coverage(mean_coverage=None,baf=None):
    '''simulate the coverage of B allele and the total coverage'''
    coverage=numpy.random.poisson(mean_coverage)
    b_allele_coverage=numpy.random.binomial(n=coverage,p=baf)
    return [coverage,b_allele_coverage]

if __name__ == '__main__':

    parse=argparse.ArgumentParser(description='Generate snvs on a coalescent tree in newick format')
    parse.add_argument('-t','--tree',required=True,help='a tree in newick format')
    default=300
    parse.add_argument('-r','--snv_rate',type=float,default=default,help='the muation rate for generating SNVs [{}]'.format(default))
    default=3
    parse.add_argument('-R','--cnv_rate',type=float,default=default,help='the muation rate for generating CNVs [{}]'.format(default))
    default=0.5
    parse.add_argument('-d','--del_prob',type=int,default=default,help='the probability to be deletion for a CNV mutation [{}]'.format(default))
    default=0.001
    parse.add_argument('-l','--cnv_length_lambda',type=float,default=default,help='the lambda of CNVs length [{}]'.format(default))
    default=0.01
    parse.add_argument('-L','--cnv_length_max',type=float,default=default,help='the maximium of CNVs length [{}]'.format(default))
    default=5
    parse.add_argument('-c','--copy_max',type=int,default=default,help='the maximium ADDITIONAL copy of a CNV [{}]'.format(default))
    default=2
    parse.add_argument('-p','--ploid',type=int,default=default,help='the ploid to simulate [{}]'.format(default))
    default=50
    parse.add_argument('-D','--depth',type=int,default=default,help='the mean depth for simulating coverage data [{}]'.format(default))
    default=None
    parse.add_argument('-s','--random_seed',type=int,help='the seed for random random number generator [{}]'.format(default))
    default='raw.cnvs'
    parse.add_argument('-V','--cnv',type=str,default=default,help='the file to save CNVs [{}]'.format(default))
    default='depth.profile'
    parse.add_argument('-P','--depth_profile',type=str,default=default,help='the file to save depth profile [{}]'.format(default))
    default='log.txt'
    parse.add_argument('-g','--log',type=str,default=default,help='the log file [{}]'.format(default))
    default=10
    parse.add_argument('--loglevel',type=int,default=default,choices=[10,20,30,40,50],help='the log file [{}]'.format(default))
    parse.add_argument('--trunk_vars',type=str,help='the trunk variants file')
    args=parse.parse_args()

    logging.basicConfig(filename=args.log, filemode='w', format='%(levelname)s: %(message)s', level=args.loglevel)
    logging.info(' Command: %s',' '.join(sys.argv))
    if args.random_seed==None:
        seed=numpy.random.randint(4294967296) #2**32: Must be convertible to 32 bit unsigned integers.
    else:
        seed=args.random_seed
    logging.info(' Random seed: %s',seed)
    numpy.random.seed(seed)
    with open(args.tree) as input:
        for line in input:
            newick=line.rstrip()
            mytree=newick2tree(newick)
#trunk vars
            trunk_snvs={}
            trunk_dels={}
            trunk_cnvs={}
            if args.trunk_vars!=None:
                trunk_snvs,trunk_dels,trunk_cnvs=trunk_vars.classify_vars(args.trunk_vars,args.ploid,mytree.leaves_number(),mytree)

            snvs_freq,cnvs,depth_profile=mytree.snvs_freq_cnvs_profile(
                                                       ploid=args.ploid,
                                                       snv_rate=args.snv_rate,
                                                       cnv_rate=args.cnv_rate,
                                                       del_prob=args.del_prob,
                                                       cnv_length_lambda=args.cnv_length_lambda,
                                                       cnv_length_max=args.cnv_length_max,
                                                       copy_max=args.copy_max,
                                                       trunk_snvs=trunk_snvs,
                                                       trunk_dels=trunk_dels,
                                                       trunk_cnvs=trunk_cnvs,
                                                       )
            cnv_file=open(args.cnv,'w')
            for cnv in cnvs:
                cnv_file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(cnv['seg'],cnv['start'],cnv['end'],cnv['copy'],cnv['leaves_count'],cnv['pre_snvs']))

            depth_profile_file=open(args.depth_profile,'w')
            for reg in depth_profile:
                depth_profile_file.write('{}\t{}\t{}\n'.format(*reg))

            for snv in snvs_freq:
                print(*snv,sep="\t")


