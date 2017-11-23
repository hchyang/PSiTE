#!/usr/bin/env python3

import re
import pickle
import numpy
import copy
import logging
import os

class Tree:
    snv_pos={}
    def __init__(self,name=None,lens=None,left=None,right=None,top=None,snvs=None,accumulated_snvs=None,cnvs=None,accumulated_cnvs=None,C='0.0.0',nodeid=None,sim=True):
        self.name=name
        self.lens=lens
        self.left=left
        self.right=right
        self.top=top #ancestor node
#it's a list of dictionary and each dictionary contains those keys {type,start,end,mutation} of each snv that occured on its top branch
        self.snvs=snvs 
        self.accumulated_snvs=accumulated_snvs #it contains pos for all snvs on the lineage leading to that node 
#it's a list of dictionary and each dictionary contains those keys {type,start,end,copy,leaves_count,pre_snvs,new_copies} of each cnv that occured on its top branch
        self.cnvs=cnvs                         
        self.accumulated_cnvs=accumulated_cnvs 
        self.C=C
        self.nodeid=nodeid
        self.sim=sim

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
            raise ShouldNotBeHereError
        return self

    #@profile
    def add_snv_cnv(self,start=None,end=None,inherent_snvs=None,inherent_cnvs=None,
                    snv_rate=None,cnv_rate=None,del_prob=None,cnv_length_beta=None,
                    cnv_length_max=None,cn_dist_cfg=None,tstv_dist_cfg=None):
        '''
        Randomly put SNVs and CNVs on a phylogenetic tree.
        For amplifications, we will build a new tree for every new copy and use this method to add SNVs/CNVs on the new tree.
        NOTE: 1. For both SNV and CNV, the position is 0 based.
              2. For each CNVs, its start is inclusive, but its end is not: [start, end).
        '''
#rescale mutation rate according the length of the sequence
        if inherent_snvs==None:
            inherent_snvs=[]
        if inherent_cnvs==None:
            inherent_cnvs=[]
        length=end-start
        logging.debug('%s with length: %s',self.nodeid,self.lens)
        logging.debug('Structure: %s',self.tree2newick())
        self.snvs=[]
        self.cnvs=[]
        self.accumulated_snvs=[]
        self.accumulated_cnvs=[]
        if self.top == None: 
#root node, may have inherent_snvs
#FIXME: inherent_snvs/truncal_snvs should also be changed to a list of dictionaries. 
            self.snvs=inherent_snvs[:]
            self.cnvs=inherent_cnvs[:]
            self.accumulated_snvs=inherent_snvs[:]
#check inherent_dels
#            self.accumulated_dels=inherent_dels[:]
            self.accumulated_cnvs=inherent_cnvs[:]
        else:
#non-root node inherits snvs/cnvs from its top nodes 
            if self.top.accumulated_snvs != None:
                self.accumulated_snvs=self.top.accumulated_snvs[:] 
            if self.top.accumulated_cnvs != None:
                self.accumulated_cnvs=self.top.accumulated_cnvs[:]
#rescale with the length
        mutation_rate=snv_rate+cnv_rate
#skip 1. self.sim==False (the leveas of which are less or equal to prune)
#     2. self.lens=None  (the root)
#     3. mutation_rate<=0
        if self.sim and self.lens!=None and mutation_rate > 0:
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
                    new_snv=1
                    if self.accumulated_cnvs:
                        for del_start,del_end in [[cnv['start'],cnv['end']] for cnv in self.accumulated_cnvs if cnv['type']=='DEL']:
                            if del_start<=pos<del_end:
                                new_snv=0
                                break
                    if new_snv==1:
                        Tree.snv_pos[pos]=1
                        snv={'type':'SNV',
                             'start':pos,
                             'end':pos+1,
                             'mutation':numpy.random.choice(tstv_dist_cfg['form'],p=tstv_dist_cfg['prob'])}
                        self.snvs.append(snv)
                        self.accumulated_snvs.append(snv)
                        logging.debug('New SNV: %s',pos)
                        logging.debug('The length of the branch new SNV locates at: %s',self.lens)
                        logging.debug('Structure: %s',self.tree2newick())
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
                    logging.debug('Previous deletions: %s',str([[cnv['start'],cnv['end']] for cnv in self.accumulated_cnvs if cnv['type']=='DEL']))
#TODO: check here very carefully.
#We need to modify new_cnvs in place. Let's sort accumulated_dels first.
#After the sorting, all deletions in accumulated_dels should be ordered and without overlapping regions.
#Without this, there will be problems.
                    self.accumulated_cnvs.sort(key=lambda cnv: cnv['start'])
                    for cnv in new_cnvs: 
                        for del_start,del_end in [[cnv['start'],cnv['end']] for cnv in self.accumulated_cnvs if cnv['type']=='DEL']:
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
                                    raise ShouldNotBeHereError
                    logging.debug('New CNVs after adjusting to previous deletions: %s',str(new_cnvs))
                    if len(new_cnvs)==0 or len(new_cnvs[0])==0:
                        continue
########################################################################################################################
                    if numpy.random.uniform()<del_prob:
#the new cnv is a deletion
                        logging.debug('New CNVs are deletions.')
                        logging.debug('%s accumulated_snvs: %s.',self.nodeid,str([x['start'] for x in self.accumulated_snvs]))
                        logging.debug('%s snvs: %s.',self.nodeid,str([x['start'] for x in self.snvs]))
                        for del_start,del_end in new_cnvs:
#output pre_snvs to self.cnvs, so it can be used to correct the count of snvs 
                            pre_snvs=[]
#We need to use a copy of self.accumulated_snvs for the 'for loop'.
#Without that, modify this list in place will cause some element bypassed.
                            for snv in self.accumulated_snvs[:]:
                                if del_start<=snv['start']<del_end:
                                    self.accumulated_snvs.remove(snv)
                                    if snv in self.snvs:
                                        self.snvs.remove(snv)
                                    else:
                                        pre_snvs.append(snv)
                            logging.debug('pre_snvs in new DELs regions: %s.',str(pre_snvs))
                            cnv={'type':'DEL',
                                 'seg':[start,end],
                                 'start':del_start,
                                 'end':del_end,
                                 'copy':-1,
                                 'leaves_count':leaves_count,
                                 'pre_snvs':pre_snvs,
                                 'new_copies':[]}
                            self.cnvs.append(cnv)
                            self.accumulated_cnvs.append(cnv)
                    else:
#the new cnv is an amplification
                        logging.debug('New CNVs are amplifications.')
                        cnv_copy=numpy.random.choice(cn_dist_cfg['copy'],p=cn_dist_cfg['prob'])
                        for amp_start,amp_end in new_cnvs:
                            amp_length=amp_end-amp_start
#collect the old snvs on cnvs. Those snvs are the snvs on the ancestor lineage leading to segment, and locate in the segment.
                            pre_snvs=[]
                            for snv in self.accumulated_snvs:
                                if amp_start<=snv['start']<amp_end:
                                    pre_snvs.append(snv)
#collect the new copies of cnvs
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
                            cnv={'type':'AMP',
                                 'seg':[start,end],
                                 'start':amp_start,
                                 'end':amp_end,
                                 'copy':cnv_copy,
                                 'leaves_count':leaves_count,
                                 'pre_snvs':pre_snvs,
                                 'new_copies':new_copies}
                            self.cnvs.append(cnv)
                            self.accumulated_cnvs.append(cnv)
        for cnv in self.cnvs:
            if cnv['copy']>0: 
                scale=(cnv['end']-cnv['start'])/(end-start)
                for segment in cnv['new_copies']:
#For each new copy of amplification, the current node is the root node. It inherent all the snvs in pre_snvs.
#But we have compared all new cnvs with accumulated_cnvs, so no pre_cnvs will affect our new copies.
                    segment.add_snv_cnv(start=cnv['start'],end=cnv['end'],inherent_snvs=cnv['pre_snvs'],
                                        snv_rate=snv_rate*scale,cnv_rate=cnv_rate*scale,del_prob=del_prob,
                                        cnv_length_beta=cnv_length_beta,cnv_length_max=cnv_length_max,
                                        cn_dist_cfg=cn_dist_cfg,tstv_dist_cfg=tstv_dist_cfg)
#only root node have inherent_snvs and inherent_cnvs
        if self.left != None:
            self.left.add_snv_cnv(start=start,end=end,inherent_snvs=[],inherent_cnvs=[],
                                  snv_rate=snv_rate,cnv_rate=cnv_rate,del_prob=del_prob,
                                  cnv_length_beta=cnv_length_beta,cnv_length_max=cnv_length_max,
                                  cn_dist_cfg=cn_dist_cfg,tstv_dist_cfg=tstv_dist_cfg)
        if self.right != None:
            self.right.add_snv_cnv(start=start,end=end,inherent_snvs=[],inherent_cnvs=[],
                                   snv_rate=snv_rate,cnv_rate=cnv_rate,del_prob=del_prob,
                                   cnv_length_beta=cnv_length_beta,cnv_length_max=cnv_length_max,
                                   cn_dist_cfg=cn_dist_cfg,tstv_dist_cfg=tstv_dist_cfg)

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

    def all_snvs_summary(self):
        '''
        It will return a dictionary of all SNVs in the tree. 
        {pos:{'mutation':,'alt_count':},...}
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
                all_alt_count[snv['start']]={'mutation':snv['mutation'],'alt_count':self.leaves_counting()}
        if self.cnvs:
            for cnv in self.cnvs:
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        all_alt_count=merge_two_all_alt_count(all_alt_count,cp.all_snvs_summary())
                else:  #deletion
                    pre_snvs_dict={}
                    for snv in cnv['pre_snvs']:
                        pre_snvs_dict[snv['start']]={'mutation':snv['mutation'],'alt_count':-self.leaves_counting()}
                    all_alt_count=merge_two_all_alt_count(all_alt_count,pre_snvs_dict)
        if self.left!=None:
            all_alt_count=merge_two_all_alt_count(all_alt_count,self.left.all_snvs_summary())
        if self.right!=None:
            all_alt_count=merge_two_all_alt_count(all_alt_count,self.right.all_snvs_summary())
        return all_alt_count

    def nodes_vars_collect(self,chroms=None):
        '''
        It will return a dictionary of all nodes in the tree. 
        {node1:{var1,var2,...},node2:{var3,var4,...},...}
        '''
        nodes_vars={}
        nodes_vars[self.nodeid]=set()
        if self.snvs:
            for snv in self.snvs:
                var='#'.join([str(x) for x in [chroms,snv['start'],snv['end'],snv['mutation']]])
                nodes_vars[self.nodeid].add(var)
        if self.cnvs:
            for cnv in self.cnvs:
                var='#'.join([str(x) for x in [chroms,cnv['start'],cnv['end']]])
                if cnv['copy']>0:
                    var+='#+'+str(cnv['copy'])
                else:
                    var+='#'+str(cnv['copy'])
                nodes_vars[self.nodeid].add(var)
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        tmp=cp.nodes_vars_collect(chroms=chroms)
#For each new copy of amplification, its self.snvs contains previous snvs.
#Those snvs do not locate on the current node, let's remove them.
#FIXME: Acutually, this is not true. If in one branch, SNV, amplificatoin and deletion occure sequentially,
#and they overlap each other, the new SNV will not captured by the current node on main tree, but
#still captured by current node on the new copy, and will be eliminated as it's in the set pre_snvs.
                        if tmp.get(self.nodeid) and cnv['pre_snvs']:
                            for snv in cnv['pre_snvs']:
                                var='#'.join([str(x) for x in [chroms,snv['start'],snv['end'],snv['mutation']]])
                                tmp[self.nodeid].discard(var)
                        nodes_vars=merge_two_dict_set(nodes_vars,tmp)
        if self.left!=None:
            nodes_vars=merge_two_dict_set(nodes_vars,self.left.nodes_vars_collect(chroms=chroms))
        if self.right!=None:
            nodes_vars=merge_two_dict_set(nodes_vars,self.right.nodes_vars_collect(chroms=chroms))
        return nodes_vars

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
    
    def attach_info(self,attr=None,info=None):
        '''
        Put the informaton of each node in the DICTIONARY (info) onto each node.
        Will set None as the default value.
        '''
        if info==None:
            info={}
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
#for a node, if node.left.leaves_count<=tips and node.right.leaves_count<=tips, prune it into a tip node
#if node.left.leaves_count<=tips and node.right.leaves_count>tips, just prune node.left into a tip node
        if self.left==None and self.right==None:
            if self.name==None:
                self.name=self.nodeid
            if self.leaves_count<=tips:
                self.sim=False
        else:
            if self.left.leaves_count<=tips and self.right.leaves_count<=tips:
                self.left=None
                self.right=None
                if self.name==None:
                    self.name=self.nodeid
                if self.leaves_count<=tips:
                    self.sim=False
            else:
                self.left.prune(tips=tips)
                self.right.prune(tips=tips)


#Let's use this method to collect all the snvs for each leaf.
    #@profile
    def genotyping(self,genotypes=None):
        '''
        Collect the genotypes on every SNV site for each leaf.
        And modify the dictionary (genotypes) directly.
        The dictionary's data structure is:
        {leaf1:{pos1:genotype,pos2:genotype...},leaf2:{pos1:genotype,pos2:genotype...},...}
        '''
        #logging.debug('snv_genotyping: %s',self.nodeid)
        if genotypes==None:
            genotypes={}
        if self.snvs:
            for leaf in self.leaves_naming():
                if leaf not in genotypes:
                    genotypes[leaf]={}
                for pos in [snv['start'] for snv in self.snvs]:
                    if pos in genotypes[leaf]:
                        genotypes[leaf][pos]+=1
                    else:
                        genotypes[leaf][pos]=1
        if self.cnvs:
            for cnv in self.cnvs:
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        cp.genotyping(genotypes)
                else:  #deletion
                    for pos in [snv['start'] for snv in cnv['pre_snvs']]:
                        for leaf in self.leaves_naming():
                            genotypes[leaf][pos]-=1
        if self.left!=None:
            self.left.genotyping(genotypes)
        if self.right!=None:
            self.right.genotyping(genotypes)

    #@profile
    def cnv_genotyping(self,genotypes=None,parental=None):
        '''
        Collect the genotypes on every CNV site for each leaf.
        '''
        #logging.debug('cnv_genotyping: %s',self.nodeid)
        if genotypes==None:
            genotypes={}
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

#######################################
#TODO
#In order to build haplotype for each node efficiently, I will
# 1. Store more information of each SNV in a dictionary, and collect all of the SNVs
#    on the lineage leading to each tip node in accumulated_SNVs
# 2. Collect all of the SNVs on the lineage leading to each tip node in accumulated_CNVs
# 3. Traverse the whole haplotype tree and add one more pair of key:value to each CNV dictionary.
#    The key is 'haplotypes', and the value is a list of dictionaries, each dictionary is:
#    {'tip_node1':[SNVs+CNVs],'tip_node2':[SNVs+CNVs],...}
#    After this operation, all CNVs in accumulated_CNVs will be changed in place.
# 4. Build a nested data structure can be used to build haplotype reference.
#    {'start':start,'end':end,'vars':{'tip_node1':[SNVs+CNVs],'tip_node2':[SNVs+CNVs],...}}
#It's not easy to build the nested data structure:
#haptype->[vars->node->CNV->haplotypes(transformed from new_copies)->haplotype]
#Each layer has 5 sub-layers.
#######################################

    def tip_node_leaves(self,tip_leaves=None):
        '''
        Return a dictionary {tip_node1:[leaf_name1,leaf_name2...],tip_node2:[...],...}.
        '''
        if tip_leaves==None:
            tip_leaves={}
        if self.left!=None:
            self.left.tip_node_leaves(tip_leaves=tip_leaves)
        if self.right!=None:
            self.right.tip_node_leaves(tip_leaves=tip_leaves)
        if self.left==None and self.right==None:
            tip_leaves[self.nodeid]=self.leaves_naming()
        return tip_leaves

    def leaf_vars(self,start=None,end=None,tip_vars=None):
        '''
        Only return the vars on the main tree level, will not trace the vars on the new copies of each CNVs.
        '''
        if tip_vars==None:
            tip_vars={}
        if self.left!=None:
            self.left.leaf_vars(start=start,end=end,tip_vars=tip_vars)
        if self.right!=None:
            self.right.leaf_vars(start=start,end=end,tip_vars=tip_vars)
        if self.left==None and self.right==None:
            if tip_vars=={}:
                logging.debug('???!!!!start: %s; end: %s',start,end)
                tip_vars['start']=start
                tip_vars['end']=end
                tip_vars['vars']={}
            else:
                logging.debug('tip_vars: %s',str(tip_vars))
            tip_vars['vars'][self.nodeid]=self.accumulated_snvs+self.accumulated_cnvs
            tip_vars['vars'][self.nodeid].sort(key=lambda var:var['start'])
        logging.debug('%s',self.nodeid)
        logging.debug('start: %s; end: %s',start,end)
        return tip_vars
    
    def add_cnvs_haps_key(self):
        '''
        In this method, I will add haplotypes to each CNV on the tree.
        This is just for CNV. For the main original chromosome, I have to add the haplotypes manually.
        I will do this in construct_leaf_haplotype.
        '''
        for cnv in self.cnvs:
            if cnv['type']=='AMP':
                cnv['haplotypes']=[]
                for copy in cnv['new_copies']:
                    cnv['haplotypes'].append(copy.leaf_vars(start=cnv['start'],end=cnv['end']))
                    copy.add_cnvs_haps_key()
        if self.left!=None:
            self.left.add_cnvs_haps_key()
        if self.right!=None:
            self.right.add_cnvs_haps_key()

    def construct_leaf_haplotype(self,start=None,end=None):
        '''
        I will collect all vars (snvs+cnvs) for each tip node here.
        This function will apply to each haplotype.
        The function will fill a dictionary with the structure:
        {'start':start,'end':end,'vars':{'tip_node1':[SNVs+CNVs],'tip_node2':[SNVs+CNVs],...}}
        !!!In this structure, each haplotype of each haplotypes in CNVs will have the same structure as above.!!!
        '''
        self.add_cnvs_haps_key()
        try:
            logging.debug('tip_vars : %s',tip_vars)
        except NameError:
            logging.debug('not exist yet')
        leaf_haplotype=self.leaf_vars(start=start,end=end)
        return leaf_haplotype

    #@profile
    def snvs_freq_cnvs_profile(self,parental=None,snv_rate=None,cnv_rate=None,del_prob=None,
                               cnv_length_beta=None,cnv_length_max=None,cn_dist_cfg=None,tstv_dist_cfg=None,
                               trunk_snvs=None,trunk_cnvs=None,purity=None,
                               length=None,chain=None,chroms=None):
        '''
        Produce the true frequency of SNVs in the samples.
        It's a warpper for generating SNVs/CNVs on a tree and summarize their frequency.
        '''
        all_cnvs=[]
        nodes_vars={}
        all_snvs_alt_counts={}
        all_snvs_alt_freq=[]
#leaf_snv_alts is a hash of hash, {leaf1:{pos1:genotype,pos2:genotype...},leaf2:{pos1:genotype,pos2:genotype...},...}
        leaf_snv_alts={}
        leaf_cnvs={}
        leaf_cnvs_pos_changes={}
        ploidy=len(parental)

        background=self.leaves_counting()*ploidy
        logging.debug('Your tree is: %s',self.tree2newick())
        hap_cnvs=[]
#collect all snvs and cnvs
        for i in range(ploidy):
            logging.info(' Simulate tree %s (total: %s)',i+1,ploidy)
            hap_tree=copy.deepcopy(self)
            hap_trunk_snvs=trunk_snvs.get(i,[])
            hap_trunk_cnvs=trunk_cnvs.get(i,[])
            hap_tree.add_snv_cnv(start=0,end=length,inherent_snvs=hap_trunk_snvs,
                inherent_cnvs=hap_trunk_cnvs,snv_rate=snv_rate,
                cnv_rate=cnv_rate,del_prob=del_prob,cnv_length_beta=cnv_length_beta,
                cnv_length_max=cnv_length_max,cn_dist_cfg=cn_dist_cfg,tstv_dist_cfg=tstv_dist_cfg)

            all_snvs_alt_counts.update(hap_tree.all_snvs_summary())

            hap_cnvs.append(hap_tree.all_cnvs_collect())
            all_cnvs.extend(hap_cnvs[-1])
            nodes_vars=merge_two_dict_set(nodes_vars,hap_tree.nodes_vars_collect(chroms=chroms))

            hap_tree.genotyping(genotypes=leaf_snv_alts)
            hap_tree.cnv_genotyping(genotypes=leaf_cnvs,parental=parental[i])
            if chain!=None:
                leaf_haplotype=hap_tree.construct_leaf_haplotype(start=0,end=length)
                logging.debug('Haplotypes: %s',leaf_haplotype)
                output_leaf_haplotype(leaf_haplotype=leaf_haplotype,directory=chain,chroms=chroms,haplotype=i,parental=parental[i])

        all_snvs_pos=sorted(all_snvs_alt_counts.keys())

#construct cnv profile list, assuming the whole region start with 0 and end with length
        all_cnvs.sort(key=lambda cnv: cnv['start'])
        cnvs_pos_changes=cnvs2pos_changes(cnvs=all_cnvs,length=length,background=background)

        hap_local_copy_for_all_snvs=hap_local_leaves(positions=all_snvs_pos,
            hap_cnvs=hap_cnvs,length=length,background=self.leaves_counting(),ploidy=ploidy)

#construct cnv profile for each hap_tree, assuming the whole region start with 0 and end with length
#translate cnvs_pos_changes to cnv_profile before its changing
        cnv_profile=pos_changes2region_profile(cnvs_pos_changes)

        region_mean_ploidy=0
        normal_dosage=background*(1-purity)/purity
        for pos in all_snvs_pos:
            while pos>=cnvs_pos_changes[0][0]:
                region_mean_ploidy+=cnvs_pos_changes.pop(0)[1]
#adjust SNVs' frequency by taking the normal cells into account 
            all_snvs_alt_freq.append([pos,all_snvs_alt_counts[pos]['mutation'],all_snvs_alt_counts[pos]['alt_count']/(normal_dosage+region_mean_ploidy)])
#TODO: what information of CNVs should I output? logR? Frequency? Adjust CNVs' frequency?

#build genotypes for each leaf
#build genotypes on the variants of each leaf? or on variants on all leaves?
        leaf_snv_refs={}
        for leaf in self.leaves_naming():
            if leaf not in leaf_cnvs:
                leaf_cnvs[leaf]=[]
            leaf_cnvs[leaf].sort(key=lambda cnv:(cnv['start'],cnv['end']))
            leaf_cnvs_pos_changes[leaf]=cnvs2pos_changes(cnvs=leaf_cnvs[leaf],length=length,background=ploidy)
            if leaf in leaf_snv_alts:
                for pos in all_snvs_pos:
                    if pos not in leaf_snv_alts[leaf]:
                        leaf_snv_alts[leaf][pos]=0
            else:
                leaf_snv_alts[leaf]={}
                for pos in all_snvs_pos:
                    leaf_snv_alts[leaf][pos]=0
            region_mean_ploidy=0
            leaf_snv_refs[leaf]={}
#            print(leaf)
#            print(leaf_cnvs_pos_changes[leaf])
            for pos in all_snvs_pos:
                while pos>=leaf_cnvs_pos_changes[leaf][0][0]:
                    region_mean_ploidy+=leaf_cnvs_pos_changes[leaf].pop(0)[1]
                leaf_snv_refs[leaf][pos]=region_mean_ploidy-leaf_snv_alts[leaf][pos]

        return all_snvs_alt_freq,all_cnvs,cnv_profile,nodes_vars,leaf_snv_alts,leaf_snv_refs,leaf_cnvs,hap_local_copy_for_all_snvs

    def tree2newick(self,lens=False,attrs=None):
        '''
        Convert tree structure to string in Newick/NHX format.
        '''
        newick_str=''

        if self.left!=None:
            newick_str+='(' + self.left.tree2newick(lens=lens,attrs=attrs)
        if self.name==None:
            newick_str+=','
        else:
            newick_str+=self.name
        if self.right!=None:
            newick_str+=self.right.tree2newick(lens=lens,attrs=attrs) + ')'
        if self.lens!=None and lens:
            newick_str+= ':' + self.lens
        if attrs!=None and lens==True:
            newick_str+='[&&NHX'
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

def merge_two_all_alt_count(dict1=None,dict2=None):
    '''
    Each dict has the structure: {pos:{'mutation':,'alt_count':},...}.
    I will generate a new dict, which contains all pos in dict1 and dict2,
    And for those pos have record in both dict1 and dict2, the vale of 
    'mutation' will be keep, and the value of 'alt_count' is the sum of value
    in dict1 and dict2.
    '''
    if dict1==None:
        dict1={}
    if dict2==None:
        dict2={}
    new_dict=dict1.copy()
    for key in dict2:
        if key in new_dict:
            new_dict[key]['alt_count']+=dict2[key]['alt_count']
        else:
            new_dict[key]=dict2[key]
    return new_dict
    
def merge_two_dict_set(dict1=None,dict2=None):
    '''
    It's similiar with dict.update, but for the key in both dicts,
    its new value equal the union of dict1[key] and dict2[key].
    '''
    if dict1==None:
        dict1={}
    if dict2==None:
        dict2={}
    new_dict=dict1.copy()
    for key in dict2:
        if key in new_dict:
            new_dict[key].update(dict2[key])
        else:
            new_dict[key]=dict2[key]
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
        yield 'node'+str(i)

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
                print('Check your newick tree! Or maybe there are something I do not know about newick!')
                raise ShouldNotBeHereError
        else:
            mytree=mytree.top
    return mytree
    
def simulate_sequence_coverage(mean_coverage=None,baf=None):
    '''
    Simulate the coverage of B allele and the total coverage
    '''
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

def output_leaf_haplotype(leaf_haplotype=None,directory=None,chroms=None,haplotype=None,parental=None):
    '''
    Output the variants of each tip node in the order of coordinate.
    '''
    for node in leaf_haplotype['vars'].keys():
        with open('{}/{}.genome.chain'.format(directory,node),'a') as cfg_file:
            cfg_file.write('>{}_Haplotype{} parental:{}\n'.format(chroms,haplotype,parental))
            retrieve_tip_vars(tip_vars=leaf_haplotype,tip=node,out_file=cfg_file,chroms=chroms)

def retrieve_tip_vars(tip_vars=None,tip=None,out_file=None,chroms=None):
    '''
    The data structure of tip_vars is:
    {'start':start,'end':end,'vars':{'tip_node1':[SNVs+CNVs],'tip_node2':[SNVs+CNVs],...}}
    In this structure, each copy of each CNV in CNVs will have the same structure as above.
    '''
    seq_seg=[]
    breakpoint=tip_vars['start']
    for var in tip_vars['vars'][tip]:
#There will be no SNV overlap with DEL.
#But AMP can overlap with SNV or DEL. As AMP will not change the breakpoint to its start,
#and all VARs are sorted by start, so the situation of breakpoint>start will only occure
#in AMP events.
        if var['type']=='SNV': #snv
            if var['start']>breakpoint:
                out_file.write(build_line(elements=[chroms,breakpoint,var['start'],'ref']))
                out_file.write(build_line(elements=[chroms,var['start'],var['end'],var['type'],var['mutation']]))
            elif var['start']==breakpoint:
                out_file.write(build_line(elements=[chroms,var['start'],var['end'],var['type'],var['mutation']]))
            else:
                raise ShouldNotBeHereError
            breakpoint=var['end']
        elif var['type']=='DEL': #deletion
            if var['start']>breakpoint:
                out_file.write(build_line(elements=[chroms,breakpoint,var['start'],'ref']))
                out_file.write(build_line(elements=[chroms,var['start'],var['end'],var['type']]))
            elif var['start']==breakpoint:
                out_file.write(build_line(elements=[chroms,var['start'],var['end'],var['type']]))
            else:
                raise ShouldNotBeHereError
            breakpoint=var['end']
        elif var['type']=='AMP': #amplification
            if var['start']>breakpoint:
                out_file.write(build_line(elements=[chroms,breakpoint,var['start'],'ref']))
                breakpoint=var['start']
            for haplotype in var['haplotypes']:
                retrieve_tip_vars(tip_vars=haplotype,tip=tip,out_file=out_file,chroms=chroms)
        else: 
            raise ShouldNotBeHereError
    if tip_vars['end']>breakpoint:
        out_file.write(build_line(elements=[chroms,breakpoint,tip_vars['end'],'ref']))
    elif tip_vars['end']==breakpoint:
        pass
    else:
        raise ShouldNotBeHereError

def build_line(elements=None):
    '''
    elements should be a list.
    I will join them with '\t' and append a '\n' at the tail.
    '''
    return '{}\n'.format('\t'.join([str(x) for x in elements]))

class TooManyMutationsError(Exception):
    pass

class ShouldNotBeHereError(Exception):
    pass
