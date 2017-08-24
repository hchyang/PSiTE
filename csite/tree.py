#!/usr/bin/env python3

import re
import pickle
import numpy
import copy
import logging

class Tree:
    snv_pos={}
    def __init__(self,name=None,lens=None,left=None,right=None,top=None,snvs=None,accumulated_snvs=None,cnvs=None,accumulated_dels=None,C='0.0.0',nodeid=None):
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
        self.nodeid=nodeid

    def add_node(self,node=None):
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

    #@profile
    def add_snv_cnv(self,start=None,end=None,inherent_snvs=[],inherent_dels=[],inherent_cnvs=[],
                    snv_rate=None,cnv_rate=None,del_prob=None,cnv_length_beta=None,cnv_length_max=None,cn_dist_cfg=None):
        '''
        Randomly put SNVs and CNVs on a phylogenetic tree.
        For amplifications, we will build a new tree for every new copy and use this method to add SNVs/CNVs on the new tree.
        NOTE: 1. For both SNV and CNV, the position is 0 based.
              2. For each CNVs, its start is inclusive, but its end is not: [start, end).
        '''
#rescale mutation rate according the length of the sequence
        length=end-start
        logging.debug('Node (%s) with length: %s',self.nodeid,self.lens)
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
#rescale with the length
        mutation_rate=snv_rate+cnv_rate
        if self.lens != None and mutation_rate > 0:
            snv_prob=snv_rate/(snv_rate+cnv_rate)
            mutation_waiting_times=waiting_times(span=self.lens,rate=mutation_rate)
            for waiting_t in mutation_waiting_times:
                pos=numpy.random.randint(start,end)
                if numpy.random.uniform()<snv_prob:
#snv
#make sure at most one SNV mutated on each position.
                    hit=0
                    while pos in Tree.snv_pos:
                        hit+=1
                        if hit>10:
                            raise TooManyMutationsError
                        pos=numpy.random.randint(start,end)
                    Tree.snv_pos[pos]=1
#TODO: Maybe I need to use a dictionary to save more information for each SNVs,
#like {pos:pos,tree_id:tree_id}. With that information, I can build haplotypes for sequence simulation.
                    self.snvs.append(pos)
                    self.accumulated_snvs.append(pos)
                    logging.debug('New SNV: %s',pos)
                    logging.debug('The length of the branch new SNV locates at: %s',self.lens)
                    logging.debug('Structure: %s',self.tree2newick())
                    if self.accumulated_dels:
                        for del_start,del_end in self.accumulated_dels:
                            if del_start<=pos<del_end:
                                self.snvs.pop()
                                self.accumulated_snvs.pop()
                                logging.debug('The new SNV locates in the previous deletion: [%s, %s]',del_start,del_end)
                                break
                else:
#cnvs
#if the new cnv overlap with accumulated_dels, compare it with the accumulated dels 
#and only keep those new regions.
                    cnv_start=pos
                    cnv_length=round(numpy.random.exponential(cnv_length_beta))
                    if cnv_length>cnv_length_max:
                        cnv_length=cnv_length_max
                    cnv_end=cnv_start+cnv_length
                    if cnv_end>end:
                        cnv_end=end
                    leaves_count=self.leaves_counting()
                    new_cnvs=[[cnv_start,cnv_end]]
                    logging.debug('New CNV: %s',str(new_cnvs))
                    logging.debug('Previous deletions: %s',str(self.accumulated_dels))
#TODO: check here very carefully.
#We need to modify new_cnvs in place. Let's sort accumulated_dels first.
#After the sorting, all deletions in accumulated_dels should be ordered non-overlapping regions.
#Without this, there will be problems.
                    self.accumulated_dels.sort(key=lambda deletion: deletion[0])
                    for cnv in new_cnvs: 
                        for del_start,del_end in self.accumulated_dels:
                            if cnv[0]<del_start:
                                if del_start<=cnv[1]<=del_end:
                                    cnv[1]=del_start
                                elif cnv[1]>del_end:
                                    new_cnvs.append([del_end,cnv[1]])
                                    cnv[1]=del_start
                                else:
#cnv_end < current_del_start, that means cnv_end is less than start points of all left dels
                                    break
                            elif del_start<=cnv[0]<=del_end:
                                if del_start<=cnv[1]<=del_end:
                                    new_cnvs.remove(cnv)
                                    break
                                elif cnv[1]>del_end:
                                    cnv[0]=del_end
                                else:
                                    print('Should not be here!3')
                    logging.debug('New CNVs after adjusting to previous deletions: %s',str(new_cnvs))
                    if len(new_cnvs)==0 or len(new_cnvs[0])==0:
                        continue
########################################################################################################################
                    if numpy.random.uniform()<del_prob:
#the new cnv is a deletion
                        logging.debug('New CNVs are deletions.')
                        logging.debug('Node (%s) accumulated_snvs: %s.',self.nodeid,str(self.accumulated_snvs))
                        logging.debug('Node (%s) snvs: %s.',self.nodeid,str(self.snvs))
                        for del_start,del_end in new_cnvs:
#output pre_snvs to self.cnvs, so it can be used to correct the count of snvs 
                            pre_snvs=[]
#We need to use a copy of self.accumulated_snvs for the 'for loop'.
#Without that, modify this list in place will cause some element bypassed.
                            for snv in self.accumulated_snvs[:]:
                                if del_start<=snv<del_end:
                                    self.accumulated_snvs.remove(snv)
                                    if snv in self.snvs:
                                        self.snvs.remove(snv)
                                    else:
                                        pre_snvs.append(snv)
                            logging.debug('pre_snvs in new DELs regions: %s.',str(pre_snvs))
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
                        cnv_copy=numpy.random.choice(cn_dist_cfg['copy'],p=cn_dist_cfg['prob'])
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
                                segment=Tree(name=self.name,lens=float(self.lens)-waiting_t,nodeid=self.nodeid)
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
                scale=(cnv['end']-cnv['start'])/(end-start)
                for segment in cnv['new_copies']:
                    segment.add_snv_cnv(start=cnv['start'],end=cnv['end'],inherent_snvs=cnv['pre_snvs'],
                                        snv_rate=snv_rate*scale,cnv_rate=cnv_rate*scale,del_prob=del_prob,
                                        cnv_length_beta=cnv_length_beta,cnv_length_max=cnv_length_max,cn_dist_cfg=cn_dist_cfg)
        if self.left != None:
            self.left.add_snv_cnv(start=start,end=end,inherent_snvs=[],
                                  snv_rate=snv_rate,cnv_rate=cnv_rate,del_prob=del_prob,
                                  cnv_length_beta=cnv_length_beta,cnv_length_max=cnv_length_max,cn_dist_cfg=cn_dist_cfg)
        if self.right != None:
            self.right.add_snv_cnv(start=start,end=end,inherent_snvs=[],
                                   snv_rate=snv_rate,cnv_rate=cnv_rate,del_prob=del_prob,
                                   cnv_length_beta=cnv_length_beta,cnv_length_max=cnv_length_max,cn_dist_cfg=cn_dist_cfg)

    def all_cnvs_collect(self):
        '''
        Return a list of all cnvs on the tree.
        '''
        all_cnvs=[]
        if self.cnvs:
            all_cnvs=self.cnvs[:]
            for cnv in self.cnvs:
                if cnv['copy']>0:
                    for cp in cnv['new_copies']:
                        all_cnvs.extend(cp.all_cnvs_collect())
        if self.left != None:
            all_cnvs.extend(self.left.all_cnvs_collect())
        if self.right != None:
            all_cnvs.extend(self.right.all_cnvs_collect())
        return all_cnvs

    def all_snvs_alt_allele_count(self):
        '''
        It will return a dictionary of all SNVs in the tree. 
        {pos:alt_allele_count,...}
        We will summary the allele count of all snvs on both main tree and all
        subtrees (new copy of cnv).
        There are three kind of snvs: 
        1) on the main tree 
        2) on the subtree (pre_snvs and new snvs)
        3) in deletions (pre_snvs)
        '''
        all_alt_count={}
        if self.snvs:
            for snv in self.snvs:
                all_alt_count[snv]=self.leaves_counting()
        if self.cnvs:
            for cnv in self.cnvs:
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        all_alt_count=merge_two_dict_count(all_alt_count,cp.all_snvs_alt_allele_count())
                else:  #deletion
                    pre_snvs_dict={}
                    for snv in cnv['pre_snvs']:
                        pre_snvs_dict[snv]=-self.leaves_counting()
                    all_alt_count=merge_two_dict_count(all_alt_count,pre_snvs_dict)
        if self.left!=None:
            all_alt_count=merge_two_dict_count(all_alt_count,self.left.all_snvs_alt_allele_count())
        if self.right!=None:
            all_alt_count=merge_two_dict_count(all_alt_count,self.right.all_snvs_alt_allele_count())
        return all_alt_count

    def nodes_snvs_collect(self):
        '''
        It will return a dictionary of all nodes in the tree. 
        {node1:{snv1,snv2,...},node2:{snv1,snv2,...},...}
        '''
        nodes_snvs={}
        if self.snvs:
            nodes_snvs[self.nodeid]=set(self.snvs)
        if self.cnvs:
            for cnv in self.cnvs:
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        tmp=cp.nodes_snvs_collect()
#For each new copy of amplification, its self.snvs contains previous snvs.
#Those snvs do not loate on the current node, let's remove them.
#FIXME: Acutually, this is not true. If in one branch, SNV, amplificatoin and deletion occure sequentially,
#and they overlap each other, the new SNV will not captured by the current node on main tree, but
#still captured by current node on the new copy, and will be eliminated as it's in the set pre_snvs.
                        if tmp.get(self.nodeid):
                            tmp[self.nodeid]=tmp[self.nodeid]-set(cnv['pre_snvs'])
                        nodes_snvs=merge_two_dict_set(nodes_snvs,tmp)
        if self.left!=None:
            nodes_snvs=merge_two_dict_set(nodes_snvs,self.left.nodes_snvs_collect())
        if self.right!=None:
            nodes_snvs=merge_two_dict_set(nodes_snvs,self.right.nodes_snvs_collect())
        return nodes_snvs

    def leaves_counting(self):
        '''
        After this method, all nodes will have the attribute of leaves count.
        '''
        if not hasattr(self,'leaves_count') or self.leaves_count == None:
            if self.left==None and self.right==None:
                self.leaves_count=1
            else:
                self.leaves_count=0
                if self.left!=None:
                    self.leaves_count+=self.left.leaves_counting()
                if self.right!=None:
                    self.leaves_count+=self.right.leaves_counting()
        return self.leaves_count

#FIXME: Should we store all the leaves' names in each node?
    def leaves_naming(self):
        '''
        After this method, ALL nodes will have the attribute leaves_names.
        '''
        if not hasattr(self,'leaves_names') or self.leaves_names == None:
            if self.left==None and self.right==None:
                self.leaves_names=[self.name]
            else:
                self.leaves_names=[]
                if self.left!=None:
                    self.leaves_names.extend(self.left.leaves_naming())
                if self.right!=None:
                    self.leaves_names.extend(self.right.leaves_naming())
        return self.leaves_names
    
    def attach_info(self,attr=None,info={}):
        '''
        Put the informaton of each node in the DICTIONARY (info) onto each node.
        Will set None as the default value.
        '''
        setattr(self,attr,info.get(self.nodeid))
        if self.left!=None:
            self.left.attach_info(attr,info)
        if self.right!=None:
            self.right.attach_info(attr,info)

    def prune(self,tips=None):
        '''
        Prune all branches with equal or less than the number of tips specified by the parameter tips.
        For a Tree object, it should run the leaves_counting() method before run this method.
        '''
        if self.leaves_count<=tips:
            logging.debug('Trim %s children of node (nodeid: %s)',self.leaves_count,self.nodeid)
            self.left=None
            self.right=None
            if self.name==None:
                self.name='i'+self.nodeid
        else:
            self.left.prune(tips=tips)
            self.right.prune(tips=tips)

    @profile
    def genotyping(self,genotypes={}):
        '''
        Collect the genotypes on every SNV site for each leaf.
        '''
        #logging.debug('snv_genotyping: %s',self.nodeid)
        if self.snvs:
            for leaf in self.leaves_naming():
                if leaf not in genotypes:
                    genotypes[leaf]={}
                for snv in self.snvs:
                    if snv in genotypes[leaf]:
                        genotypes[leaf][snv]+=1
                    else:
                        genotypes[leaf][snv]=1
        if self.cnvs:
            for cnv in self.cnvs:
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        cp.genotyping(genotypes)
                else:  #deletion
                    for snv in cnv['pre_snvs']:
                        for leaf in self.leaves_naming():
                            genotypes[leaf][snv]-=1
        if self.left!=None:
            self.left.genotyping(genotypes)
        if self.right!=None:
            self.right.genotyping(genotypes)

    @profile
    def cnv_genotyping(self,genotypes={},parental=None):
        '''
        Collect the genotypes on every CNV site for each leaf.
        '''
        #logging.debug('cnv_genotyping: %s',self.nodeid)
        if self.cnvs:
            for cnv in self.cnvs:
                for leaf in self.leaves_naming():
                    if leaf not in genotypes:
                        genotypes[leaf]=[]
                    else:
#In order to use cnvs2pos_changes later, we set leaves_count as 1 here.
#Actully, it makes sense, as each leaf's leaf count is 1.
                        genotypes[leaf].append({'start':cnv['start'],'end':cnv['end'],'copy':cnv['copy'],'leaves_count':1,'parental':parental})
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        cp.cnv_genotyping(genotypes=genotypes,parental=parental)
        if self.left!=None:
            self.left.cnv_genotyping(genotypes=genotypes,parental=parental)
        if self.right!=None:
            self.right.cnv_genotyping(genotypes=genotypes,parental=parental)

    @profile
    def snvs_freq_cnvs_profile(self,ploidy=None,snv_rate=None,cnv_rate=None,del_prob=None,
                               cnv_length_beta=None,cnv_length_max=None,cn_dist_cfg=None,
                               trunk_snvs=None,trunk_dels=None,trunk_cnvs=None,purity=None,length=None):
        '''
        Produce the true frequency of SNVs in the samples.
        It's a warpper for generating SNVs/CNVs on a tree and summarize their frequency.
        '''
        all_cnvs=[]
        nodes_snvs={}
        all_snvs_alt_counts=[]
        all_snvs_alt_freq=[]
#leaf_snv_alts is a hash of hash, {leaf1:{pos1:genotype,pos2:genotype...},leaf2:{pos1:genotype,pos2:genotype...},...}
        leaf_snv_alts={}
        leaf_cnvs={}
        leaf_cnvs_pos_changes={}

        background=self.leaves_counting()*ploidy
        logging.debug('Your tree is: %s',self.tree2newick())
        hap_cnvs=[]
#collect all snvs and cnvs
        for i in range(ploidy):
            logging.info(' Simulate tree %s (total: %s)',i+1,ploidy)
            hap_tree=copy.deepcopy(self)
            hap_trunk_snvs=trunk_snvs.get(i,[])
            hap_trunk_dels=trunk_dels.get(i,[])
            hap_trunk_cnvs=trunk_cnvs.get(i,[])
            hap_tree.add_snv_cnv(start=0,end=length,inherent_snvs=hap_trunk_snvs,inherent_dels=hap_trunk_dels,inherent_cnvs=hap_trunk_cnvs,
                                 snv_rate=snv_rate,cnv_rate=cnv_rate,del_prob=del_prob,
                                 cnv_length_beta=cnv_length_beta,cnv_length_max=cnv_length_max,cn_dist_cfg=cn_dist_cfg)

            all_snvs_alt_counts_dict=hap_tree.all_snvs_alt_allele_count()
            all_snvs_alt_counts.extend([[snv,all_snvs_alt_counts_dict[snv]] for snv in all_snvs_alt_counts_dict])

            hap_cnvs.append(hap_tree.all_cnvs_collect())
            all_cnvs.extend(hap_cnvs[-1])
            nodes_snvs=merge_two_dict_set(nodes_snvs,hap_tree.nodes_snvs_collect())

            hap_tree.genotyping(genotypes=leaf_snv_alts)
            hap_tree.cnv_genotyping(genotypes=leaf_cnvs,parental=i)

#construct a tree with all snvs
#FIXME: right now, it does not consider the deletion effect on pre_snvs.
        tree_with_snvs=copy.deepcopy(self)
        tree_with_snvs.attach_info(attr='new_snvs',info=nodes_snvs)
        all_snvs_alt_counts.sort(key=lambda snv: snv[0])

#construct cnv profile list, assuming the whole region start with 0 and end with length
        all_cnvs.sort(key=lambda cnv: cnv['start'])
        cnvs_pos_changes=cnvs2pos_changes(cnvs=all_cnvs,length=length,background=background)

#construct cnv profile for each hap_tree, assuming the whole region start with 0 and end with length
        hap_local_copy_for_all_snvs=hap_local_leaves(positions=[snp[0] for snp in all_snvs_alt_counts],hap_cnvs=hap_cnvs,length=length,background=self.leaves_counting(),ploidy=ploidy)

#translate cnvs_pos_changes to cnv_profile before its changing
        cnv_profile=pos_changes2region_profile(cnvs_pos_changes)

        region_mean_ploidy=0
        normal_dosage=background*(1-purity)/purity
#        print(purity)
#        print(normal_dosage)
        for snv in all_snvs_alt_counts:
            while snv[0]>=cnvs_pos_changes[0][0]:
                region_mean_ploidy+=cnvs_pos_changes.pop(0)[1]
#adjust SNVs' frequency by taking the normal cells into account 
            all_snvs_alt_freq.append([snv[0],snv[1]/(normal_dosage+region_mean_ploidy)])
#TODO: what information of CNVs should I output? logR? Frequency? Adjust CNVs' frequency?

#build genotypes for each leaf
#build genotypes on the variants of each leaf? or on variants on all leaves?
        leaf_snv_refs={}
        all_snvs_pos=[x[0] for x in all_snvs_alt_counts]
        for leaf in self.leaves_naming():
            if leaf not in leaf_cnvs:
                leaf_cnvs[leaf]=[]
            leaf_cnvs[leaf].sort(key=lambda cnv:(cnv['start'],cnv['end']))
            leaf_cnvs_pos_changes[leaf]=cnvs2pos_changes(cnvs=leaf_cnvs[leaf],length=length,background=ploidy)
            if leaf in leaf_snv_alts:
                for snv in all_snvs_pos:
                    if snv not in leaf_snv_alts[leaf]:
                        leaf_snv_alts[leaf][snv]=0
            else:
                leaf_snv_alts[leaf]={}
                for snv in all_snvs_pos:
                    leaf_snv_alts[leaf][snv]=0
            region_mean_ploidy=0
            leaf_snv_refs[leaf]={}
#            print(leaf)
#            print(leaf_cnvs_pos_changes[leaf])
            for snv in all_snvs_pos:
                while snv>=leaf_cnvs_pos_changes[leaf][0][0]:
                    region_mean_ploidy+=leaf_cnvs_pos_changes[leaf].pop(0)[1]
                leaf_snv_refs[leaf][snv]=region_mean_ploidy-leaf_snv_alts[leaf][snv]

        return all_snvs_alt_freq,all_cnvs,cnv_profile,nodes_snvs,tree_with_snvs,leaf_snv_alts,leaf_snv_refs,leaf_cnvs,hap_local_copy_for_all_snvs

    def tree2newick(self,lens=False,attrs=None):
        '''
        Convert tree structure to string in Newick/NHX format.
        '''
        newick_str=''

        if self.left != None:
            newick_str+='(' + self.left.tree2newick(lens=lens,attrs=attrs)
        if self.name==None:
            newick_str+=','
        else:
            newick_str+=self.name
        if self.right != None:
            newick_str+=self.right.tree2newick(lens=lens,attrs=attrs) + ')'
        if self.lens!=None and lens:
            newick_str+= ':' + self.lens
        if attrs!=None and lens==True:
            newick_str+='[&&NHX:W=1.0'
            for attribute in attrs:
                if getattr(self, attribute):
                    newick_str+=':{}={}'.format(attribute,getattr(self, attribute))
            newick_str+=']'
        return newick_str

    def highlight_snvs(self,snvs=None):
        '''
        Set the color of the nodes with certain SNVs as 255.0.0.
        '''
        if self.left != None:
            self.left.highlight_snvs(snvs)
        if self.right != None:
            self.right.highlight_snvs(snvs)
        if self.new_snvs !=None and self.new_snvs.intersection(snvs):
            self.C='255.0.0'
            logging.debug('Highlight node %s (with leaves %s) because of variants:\n%s',
                self.nodeid,self.leaves_counting(),self.new_snvs.intersection(snvs))

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

def merge_two_dict_count(dict1={},dict2={}):
    '''
    It's similiar with dict.update, but for the key in both dicts,
    its new value equal the sum of dict1[key] and dict2[key].
    '''
    new_dict={}
    for key in set.union(set(dict1),set(dict2)):
        new_dict[key]=0
        if key in dict1:
            new_dict[key]+=dict1[key]
        if key in dict2:
            new_dict[key]+=dict2[key]
    return new_dict
    
def merge_two_dict_set(dict1={},dict2={}):
    '''
    It's similiar with dict.update, but for the key in both dicts,
    its new value equal the union of dict1[key] and dict2[key].
    '''
    new_dict={}
    for key in set.union(set(dict1),set(dict2)):
        new_dict[key]=set()
        if key in dict1:
            new_dict[key]=new_dict[key].union(dict1[key])
        if key in dict2:
            new_dict[key]=new_dict[key].union(dict2[key])
    return new_dict
    
def cnvs2pos_changes(cnvs=None,length=None,background=None):
    '''
    Return a list of lists. Each sublist contain two elements. The first is the postion, and the second 
    is the copy number CHANGES across all the samples between that positon and the next position.
    [[pos,relative_copy_number_change],...]
    '''
    pos_changes=[[0,background],[length,-background]]
    for cnv in cnvs:
#        print(cnv)
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

def node_id():
    i=0
    while True:
        i+=1
        yield str(i)

#@profile
def newick2tree(newick=None):
    leaf_name_re=re.compile('^\w+:')
    lens_re=re.compile('^:[0-9.]+')
    brushwood=newick.split(',')
    node_id_gen=node_id()
    for branch in brushwood:
        while branch != '':
            if lens_re.match(branch):
                m=lens_re.match(branch)
                index=m.span()
                lens=branch[1:index[1]]
                branch=branch[index[1]:]
                mytree.lens=lens
            elif branch.startswith(')'):
                branch=branch[1:]
                mytree=mytree.top
            elif branch.startswith('('):
                if 'mytree' in vars():
                    mytree=mytree.add_node(Tree(nodeid=node_id_gen.__next__()))
                else:
                    mytree=Tree(nodeid=node_id_gen.__next__())
                branch=branch[1:]
            elif leaf_name_re.match(branch):
                m=leaf_name_re.match(branch)
                index=m.span()
                leaf_name=branch[:index[1]-1]
                branch=branch[index[1]-1:]
                mytree=mytree.add_node(Tree(name=leaf_name,nodeid=node_id_gen.__next__()))
            elif branch==';':
                break
            else:
                print('Should not be here!4')
                print('Check your newick tree! Or maybe there are something I do not know about newick!')
        else:
            mytree=mytree.top
    return mytree
    
def simulate_sequence_coverage(mean_coverage=None,baf=None):
    '''simulate the coverage of B allele and the total coverage'''
    coverage=numpy.random.poisson(mean_coverage)
    b_allele_coverage=numpy.random.binomial(n=coverage,p=baf)
    return [coverage,b_allele_coverage]

def hap_local_leaves(positions=None,hap_cnvs=None,length=None,background=None,ploidy=None):
    '''
    Calculates the local copy on each haplotype for each snvs.
    Return a list of lists, each sublist is have these elements:
    Position
    local_copy_on_hap1
    local_copy_on_hap2
    ...
    local_copy_on_hapN
    '''
    all_pos_local_copy=[]
    hap_cnvs_pos_changes=[]
    hap_local_copy=[]
    for i in range(ploidy):
        hap_cnvs[i].sort(key=lambda cnv: cnv['start'])
        hap_cnvs_pos_changes.append(cnvs2pos_changes(cnvs=hap_cnvs[i],length=length,background=background))
        hap_local_copy.append(0)
    for pos in positions:
        all_pos_local_copy.append([pos])
        for i in range(ploidy):
            while pos>=hap_cnvs_pos_changes[i][0][0]:
                hap_local_copy[i]+=hap_cnvs_pos_changes[i].pop(0)[1]
            all_pos_local_copy[-1].append(hap_local_copy[i])
    return all_pos_local_copy

class TooManyMutationsError(Exception):
    pass
