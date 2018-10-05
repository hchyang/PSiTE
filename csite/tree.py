#!/usr/bin/env python3

import re
import pickle
import numpy
import copy
import logging
import os

class Tree:
    snv_pos=set()
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
        logging.debug('Structure: %s',self.tree2nhx())
        self.snvs=[]
        self.cnvs=[]
        self.accumulated_snvs=[]
        self.accumulated_cnvs=[]
        if self.top == None: 
#root node, may have inherent_snvs
            self.snvs=inherent_snvs[:]
            self.cnvs=inherent_cnvs[:]
            self.accumulated_snvs=inherent_snvs[:]
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
                            raise TooManyHitsOnSamePosError
                        pos=numpy.random.randint(start,end)
                    new_snv=True
                    if self.accumulated_cnvs:
                        for del_start,del_end in [[cnv['start'],cnv['end']] for cnv in self.accumulated_cnvs if cnv['type']=='DEL']:
                            if del_start<=pos<del_end:
                                new_snv=False
                                break
                    if new_snv:
                        Tree.snv_pos.add(pos)
                        snv={'type':'SNV',
                             'start':pos,
                             'end':pos+1,
                             'mutation':numpy.random.choice(tstv_dist_cfg['form'],p=tstv_dist_cfg['prob'])}
                        self.snvs.append(snv)
                        self.accumulated_snvs.append(snv)
                        logging.debug('New SNV: %s',pos)
                        logging.debug('The length of the branch new SNV locates at: %s',self.lens)
                        logging.debug('Structure: %s',self.tree2nhx())
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
#We need to modify new_cnvs in place. Let's sort cnvs in accumulated_cnvs first.
#After the sorting, all deletions in accumulated_cnvs should be ordered and without overlapping regions.
#Without this, there will be problems.
                    self.accumulated_cnvs.sort(key=lambda cnv: (cnv['start'],cnv['end']))
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
                                 'pre_snvs':{0:pre_snvs},
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
                                segment=Tree(name=self.name,lens=self.lens-waiting_t,nodeid=self.nodeid,sim=self.sim)
                                if hasattr(self,'sectors'):
                                    segment.sectors=self.sectors
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
                                 'pre_snvs':{},
                                 'new_copies':new_copies}
                            for i in range(cnv_copy): cnv['pre_snvs'][i+1]=pre_snvs
                            self.cnvs.append(cnv)
                            self.accumulated_cnvs.append(cnv)
        for cnv in self.cnvs:
            if cnv['copy']>0: 
                scale=(cnv['end']-cnv['start'])/(end-start)
                for i in range(len(cnv['new_copies'])):
                    segment=cnv['new_copies'][i]
#For each new copy of amplification, the current node is the root node. It inherent all the snvs in pre_snvs.
#But we have compared all new cnvs with accumulated_cnvs, so no pre_cnvs will affect our new copies.
                    segment.add_snv_cnv(start=cnv['start'],end=cnv['end'],inherent_snvs=cnv['pre_snvs'][i+1],
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

    def all_cnvs_collect(self,sector=None):
        '''
        Return a list of all cnvs on the tree.
        '''
        all_cnvs=[]
        if self.cnvs:
            all_cnvs=[]
            for cnv in self.cnvs:
                cnv_cp=cnv.copy()
                cnv_cp['node']=self.nodeid
                cnv_cp['leaves_count']=self.sectors[sector]
                all_cnvs.append(cnv_cp)
            for cnv in self.cnvs:
                if cnv['copy']>0:
                    for cp in cnv['new_copies']:
                        all_cnvs.extend(cp.all_cnvs_collect(sector=sector))
        if self.left != None:
            all_cnvs.extend(self.left.all_cnvs_collect(sector=sector))
        if self.right != None:
            all_cnvs.extend(self.right.all_cnvs_collect(sector=sector))
        return all_cnvs

    def all_snvs_summary(self,sector=None):
        '''
        It will return a dictionary of all SNVs on the tree. 
        {pos:{'mutation':xxx,'alt_count':xxx,'node':xxx},...}
        We will summary the allele count of all snvs on the main tree and all
        subtrees (new copy of cnv).
        There are three kinds of snvs: 
        1) on the main tree 
        2) on the subtree (pre_snvs and new snvs)
        3) in deletions (pre_snvs)
        '''
        all_alt_count={}
        if self.snvs:
            for snv in self.snvs:
                all_alt_count[snv['start']]={'mutation':snv['mutation'],'alt_count':self.sectors[sector],'node':self.nodeid}
        if self.cnvs:
            for cnv in self.cnvs:
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        all_alt_count=merge_two_all_alt_count(all_alt_count,cp.all_snvs_summary(sector=sector))
                else:  #deletion
                    pre_snvs_dict={}
                    for snv in cnv['pre_snvs'][0]:
                        pre_snvs_dict[snv['start']]={'mutation':snv['mutation'],'alt_count':-self.sectors[sector],'node':self.nodeid}
                    all_alt_count=merge_two_all_alt_count(all_alt_count,pre_snvs_dict)
        if self.left!=None:
            all_alt_count=merge_two_all_alt_count(all_alt_count,self.left.all_snvs_summary(sector=sector))
        if self.right!=None:
            all_alt_count=merge_two_all_alt_count(all_alt_count,self.right.all_snvs_summary(sector=sector))
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
                    for i in range(len(cnv['new_copies'])):
                        cp=cnv['new_copies'][i]
                        tmp=cp.nodes_vars_collect(chroms=chroms)
#For each new copy of amplification, its self.snvs contains previous snvs.
#Those snvs do not locate on the current node, let's remove them.
                        if tmp.get(self.nodeid) and cnv['pre_snvs'][i+1]:
                            for snv in cnv['pre_snvs'][i+1]:
                                if self.top and snv in self.top.accumulated_snvs:
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
    
    def collect_tipnodes(self):
        '''
        After this method, ALL nodes will have the attribute tipnodes.
        '''
        if not hasattr(self,'tipnodes') or self.tipnodes == None:
            if self.left==None and self.right==None:
                self.tipnodes=[self.nodeid]
            else:
                self.tipnodes=[]
                if self.left!=None:
                    self.tipnodes.extend(self.left.collect_tipnodes())
                if self.right!=None:
                    self.tipnodes.extend(self.right.collect_tipnodes())
        return self.tipnodes
    
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

    def collect_leaves_and_trim(self,tipnode_leaves=None,sectors=None):
        '''
        NOTE: For a Tree object, you should run the leaves_counting() method on it before running this method.
        1) Prune all branches with LESS than the number of leaves specified by the prune_n in each sector.
        For a node, if node.left.leaves_count<cutoff and node.right.leaves_count<cutoff, prune it into a tip node.
        If node.left.leaves_count<cutoff and node.right.leaves_count>=cutoff, just prune node.left into a tip node,
        and then check the node.right's left and right nodes.
        2) Mark a tipnode with sim=False, if its leaves_count is less than the cutoffs of ALL sectors. 
        3) In order to uniform the names in --nhx output, we rename each tipnode with its nodeid.
        4) Add a dictionary to each node to save the sector information. The structure of this dictionary is:
        {sector1:number_of_sector1_cells_in_all_leaves_of_this_node,
         sector2:number_of_sector2_cells_in_all_leaves_of_this_node,
         ...
        }
        '''
        if not hasattr(self,'sectors'):
            self.sectors={}
        self.sim=False
        if self.left==None and self.right==None:
            tipnode_leaves[self.nodeid]=self.leaves_naming()
            self.name=self.nodeid
            for sector in sectors:
                cells=sectors[sector]['members']
                cutoff=sectors[sector]['prune_n']
                focal_cells=cells.intersection(self.leaves_names)
                self.sectors[sector]=len(focal_cells)
                if self.sim==False and len(focal_cells)>=cutoff:
                    self.sim=True
        else:
            for sector in sectors:
                cells=sectors[sector]['members']
                cutoff=sectors[sector]['prune_n']
                focal_cells=cells.intersection(self.leaves_naming())
                self.sectors[sector]=len(focal_cells)
                if self.sim==False and (len(cells.intersection(self.left.leaves_names))>=cutoff or len(cells.intersection(self.right.leaves_names))>=cutoff):
                    self.left.collect_leaves_and_trim(tipnode_leaves=tipnode_leaves,sectors=sectors)
                    self.right.collect_leaves_and_trim(tipnode_leaves=tipnode_leaves,sectors=sectors)
                    self.sim=True
            if self.sim==False:
                self.left=None
                self.right=None
                tipnode_leaves[self.nodeid]=self.leaves_naming()
                self.name=self.nodeid
                for sector in sectors:
                    cells=sectors[sector]['members']
                    cutoff=sectors[sector]['prune_n']
                    focal_cells=cells.intersection(self.leaves_names)
                    self.sectors[sector]=len(focal_cells)
                    if self.sim==False and len(focal_cells)>=cutoff:
                        self.sim=True

    def nodes_ccf(self,sectors_size=None,nodes_ccf=None):
        '''
        NOTE: For a Tree object, you should run the prune() method on it before running this method.
        Calculate the CCF (cancer cell fraction) of each node in each sector.
        '''
        nodes_ccf[self.nodeid]={}
        for sector in sectors_size:
            ncells=self.sectors[sector]
            nodes_ccf[self.nodeid][sector]=ncells/sectors_size[sector]
        if self.left!=None:
            self.left.nodes_ccf(sectors_size=sectors_size,nodes_ccf=nodes_ccf)
        if self.right!=None:
            self.right.nodes_ccf(sectors_size=sectors_size,nodes_ccf=nodes_ccf)

    def collect_sectors_nodes(self,sectors=None):
        '''
        Collect the nodes that are visible to each sector.
        After this method, the sectors will have a item 'nodes':{node1,node2,...}.
        '''
        for sector in sectors:
            cells=sectors[sector]['members']
            if cells.intersection(self.leaves_naming()):
                try:
                    sectors[sector]['nodes'].add(self.nodeid)
                except KeyError:
                    sectors[sector]['nodes']={self.nodeid}
        if self.left!=None:
            self.left.collect_sectors_nodes(sectors=sectors)
        if self.right!=None:
            self.right.collect_sectors_nodes(sectors=sectors)

    def prune(self,sectors=None):
        '''
        After this method, the root node will have an attribute tipnode_leaves,
        which is a dictionary in the form of {tipnode1:[leaf1,leaf2,...],tipnode2:[leaf3,...],...}
        '''
        if not hasattr(self,'tipnode_leaves'):
            tipnode_leaves={}
            self.collect_leaves_and_trim(tipnode_leaves=tipnode_leaves,sectors=sectors)
            self.tipnode_leaves=tipnode_leaves
        else:
            raise TreePruneError('Can not prune a tree which is pruned before!')

#Let's use this method to collect all the snvs for each tipnode.
    #@profile
    def genotyping(self,genotypes=None):
        '''
        Collect the genotypes on every SNV site for each tipnode.
        And modify the dictionary (genotypes) directly.
        The dictionary's data structure is:
        {tipnode1:{pos1:genotype,pos2:genotype...},tipnode2:{pos1:genotype,pos2:genotype...},...}
        The genotype here means the count of alternative alleles of each SNV.
        '''
        #logging.debug('snv_genotyping: %s',self.nodeid)
        if genotypes==None:
            genotypes={}
        if self.snvs:
            for tipnode in self.collect_tipnodes():
                if tipnode not in genotypes:
                    genotypes[tipnode]={}
                for pos in [snv['start'] for snv in self.snvs]:
                    if pos in genotypes[tipnode]:
                        genotypes[tipnode][pos]+=1
                    else:
                        genotypes[tipnode][pos]=1
        if self.cnvs:
            for cnv in self.cnvs:
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        cp.genotyping(genotypes)
                else:  #deletion
                    for pos in [snv['start'] for snv in cnv['pre_snvs'][0]]:
                        for tipnode in self.collect_tipnodes():
                            genotypes[tipnode][pos]-=1
        if self.left!=None:
            self.left.genotyping(genotypes)
        if self.right!=None:
            self.right.genotyping(genotypes)

    #@profile
    def cnv_genotyping(self,genotypes=None,parental=None):
        '''
        Collect the genotypes on every CNV site for each tipnode.
        '''
        #logging.debug('cnv_genotyping: %s',self.nodeid)
        if genotypes==None:
            genotypes={}
        if self.cnvs:
            for cnv in self.cnvs:
                for tipnode in self.collect_tipnodes():
                    if tipnode not in genotypes:
                        genotypes[tipnode]=[]
#set leaves_count=1 here, as a tipnode is a representative of each one of the leaves under it
                    genotypes[tipnode].append({'start':cnv['start'],'end':cnv['end'],'copy':cnv['copy'],'leaves_count':1,'parental':parental})
                if cnv['copy']>0: #amplification
                    for cp in cnv['new_copies']:
                        cp.cnv_genotyping(genotypes=genotypes,parental=parental)
        if self.left!=None:
            self.left.cnv_genotyping(genotypes=genotypes,parental=parental)
        if self.right!=None:
            self.right.cnv_genotyping(genotypes=genotypes,parental=parental)

#######################################
#In order to build haplotype for each tipnode efficiently, I will
# 1. Store more information of each SNV in a dictionary, and collect all of the SNVs
#    on the lineage leading to each tipnode in accumulated_snvs
# 2. Collect all of the SNVs on the lineage leading to each tipnode in accumulated_CNVs
# 3. Traverse the whole haplotype tree and add one more pair of key:value to each CNV dictionary.
#    The key is 'haplotypes', and the value is a list of dictionaries, each dictionary is:
#    {'tip_node1':[SNVs+CNVs],'tip_node2':[SNVs+CNVs],...}
#    After this operation, all CNVs in accumulated_CNVs will be changed in place.
# 4. Build a nested data structure can be used to build haplotype reference.
#    {'start':start,'end':end,'vars':{'tip_node1':[SNVs+CNVs],'tip_node2':[SNVs+CNVs],...}}
#It's not easy to build the nested data structure:
#haptype->[vars->tipnode->CNV->haplotypes(transformed from new_copies)->haplotype]
#Each layer has 5 sub-layers.
#######################################

    def tipnode_accumulated_vars(self,start=None,end=None,tip_vars=None):
        '''
        Only return the vars on the main tree level, will not trace the vars on the new copies of each CNVs.
        '''
        if tip_vars==None:
            tip_vars={}
        if self.left!=None:
            self.left.tipnode_accumulated_vars(start=start,end=end,tip_vars=tip_vars)
        if self.right!=None:
            self.right.tipnode_accumulated_vars(start=start,end=end,tip_vars=tip_vars)
        if self.left==None and self.right==None:
            if tip_vars=={}:
                tip_vars['start']=start
                tip_vars['end']=end
                tip_vars['vars']={}
            tip_vars['vars'][self.nodeid]=self.accumulated_snvs+self.accumulated_cnvs
            tip_vars['vars'][self.nodeid].sort(key=lambda var:var['start'])
        return tip_vars
    
    def add_haps2cnv(self):
        '''
        In this method, I will add haplotypes to each CNV on the tree.
        '''
        for cnv in self.cnvs:
            if cnv['type']=='AMP':
                cnv['haplotypes']=[]
                for copy in cnv['new_copies']:
                    cnv['haplotypes'].append(copy.tipnode_accumulated_vars(start=cnv['start'],end=cnv['end']))
                    copy.add_haps2cnv()
        if self.left!=None:
            self.left.add_haps2cnv()
        if self.right!=None:
            self.right.add_haps2cnv()

    def construct_tipnode_hap(self,start=None,end=None):
        '''
        I will collect all vars (snvs+cnvs) for each tipnode here.
        This function will apply to each haplotype.
        The function will fill a dictionary with the structure:
        {'start':start,'end':end,'vars':{'tip_node1':[SNVs+CNVs],'tip_node2':[SNVs+CNVs],...}}
        !!!In this structure, each haplotype of each haplotypes in CNVs will have the same structure as above.!!!
        '''
        self.add_haps2cnv()
        tipnode_hap=self.tipnode_accumulated_vars(start=start,end=end)
        return tipnode_hap

    #@profile
    def snvs_freq_cnvs_profile(self,parental=None,snv_rate=None,cnv_rate=None,del_prob=None,
                               cnv_length_beta=None,cnv_length_max=None,cn_dist_cfg=None,tstv_dist_cfg=None,
                               trunk_snvs=None,trunk_cnvs=None,length=None,
                               normal_dosage=None,chain=None,chroms=None,sectors=None,wholeT=None):
        '''
        Produce the true frequency of SNVs in the samples.
        It's a warpper for generating SNVs/CNVs on a tree and summarize their frequency.
        '''
        all_cnvs={}
        nodes_vars={}
        all_snvs_alt_counts={}
        all_snvs_alt_freq=[]
#tipnode_snv_alts is a hash of hash, {tipnode1:{pos1:genotype,pos2:genotype...},tipnode2:{pos1:genotype,pos2:genotype...},...}
        tipnode_snv_alts={}
        tipnode_cnvs={}
        tipnode_cnvs_pos_changes={}
        ploidy=len(parental)

        background=self.leaves_counting()*ploidy
        logging.debug('Your tree is: %s',self.tree2nhx())
#I used haps_cnvs to calculate the count number of each parental copies before.
#We do not need this feature anymore.
#        haps_cnvs=[]
#In order to avoid two SNVs occuring at the same position, I stored all the SNVs
#of each chromosome (multiple haplotype) to the set Tree.snv_pos.
        Tree.snv_pos=set()
        for snvs in trunk_snvs.values():
            Tree.snv_pos.update([snv['start'] for snv in snvs])
#collect all snvs and cnvs
        for i in range(ploidy):
            logging.info(' Simulate haplotype %s (total: %s)',i+1,ploidy)
            hap_tree=copy.deepcopy(self)
            hap_trunk_snvs=trunk_snvs.get(i,[])
            hap_trunk_cnvs=trunk_cnvs.get(i,[])
            hap_tree.add_snv_cnv(start=0,end=length,inherent_snvs=hap_trunk_snvs,
                inherent_cnvs=hap_trunk_cnvs,snv_rate=snv_rate,
                cnv_rate=cnv_rate,del_prob=del_prob,cnv_length_beta=cnv_length_beta,
                cnv_length_max=cnv_length_max,cn_dist_cfg=cn_dist_cfg,tstv_dist_cfg=tstv_dist_cfg)

#Update the dictionary all_snvs_alt_counts here
#There will not be two snps occure on the same position of different haplotype,
#except they are specified by users in trunk_vars 
            hap_trunk_snvs_pos=[snv['start'] for snv in hap_trunk_snvs]
            for sector in sectors.keys():
                if sector not in all_snvs_alt_counts:
                    all_snvs_alt_counts[sector]={}
                for pos,info in hap_tree.all_snvs_summary(sector=sector).items():
                    if pos in all_snvs_alt_counts[sector]:
                        if pos in hap_trunk_snvs_pos:
                            all_snvs_alt_counts[sector][pos]['alt_count']+=info['alt_count']
                        else:
                            raise ShouldNotBeHereError
                    else:
                        all_snvs_alt_counts[sector][pos]=info

                if sector not in all_cnvs:
                    all_cnvs[sector]=[]
                haplotype_cnvs=hap_tree.all_cnvs_collect(sector=sector)
                all_cnvs[sector].extend(haplotype_cnvs)

            nodes_vars=merge_two_dict_set(nodes_vars,hap_tree.nodes_vars_collect(chroms=chroms))
            hap_tree.genotyping(genotypes=tipnode_snv_alts)
            hap_tree.cnv_genotyping(genotypes=tipnode_cnvs,parental=parental[i])
            if chain!=None:
                tipnode_hap=hap_tree.construct_tipnode_hap(start=0,end=length)
                logging.debug('Haplotypes: %s',tipnode_hap)
                output_tipnode_hap(tipnode_hap=tipnode_hap,directory=chain,chroms=chroms,haplotype=i,parental=parental[i])

        all_snvs_pos=sorted(all_snvs_alt_counts[wholeT].keys())

#output true freq for multi-sectoring data
        for sector,info in sectors.items():
            sector_cnvs=all_cnvs[sector]
            sector_snvs_alt_counts=all_snvs_alt_counts[sector]

            sector_snvs_pos=sorted(sector_snvs_alt_counts.keys())

            sector_cnvs_pos_changes=cnvs2pos_changes(cnvs=sector_cnvs,length=length,background=len(info['members'])*ploidy)
            sector_cnv_profile=pos_changes2region_profile(sector_cnvs_pos_changes)
            sector_normal_dosage=info['normal_dosage']
            sector_local_tumor_dosage=0
            sector_snvs_alt_freq=[]
            for pos in sector_snvs_pos:
                while pos>=sector_cnvs_pos_changes[0][0]:
                    sector_local_tumor_dosage+=sector_cnvs_pos_changes.pop(0)[1]
                sector_snvs_alt_freq.append([pos,
                    sector_snvs_alt_counts[pos]['mutation'],
                    sector_snvs_alt_counts[pos]['alt_count']/(sector_normal_dosage+sector_local_tumor_dosage)])
            info['snvs_alt_freq']=sector_snvs_alt_freq
            info['cnvs']=sector_cnvs
            info['cnv_profile']=sector_cnv_profile

#calculate the number of reference alleles of each SNV for each tipnode
        tipnode_snv_refs={}
        for tipnode in self.collect_tipnodes():
            if tipnode not in tipnode_cnvs:
                tipnode_cnvs[tipnode]=[]
            tipnode_cnvs[tipnode].sort(key=lambda cnv:(cnv['start'],cnv['end']))
            tipnode_cnvs_pos_changes[tipnode]=cnvs2pos_changes(cnvs=tipnode_cnvs[tipnode],length=length,background=ploidy)
            if tipnode in tipnode_snv_alts:
                for pos in all_snvs_pos:
                    if pos not in tipnode_snv_alts[tipnode]:
                        tipnode_snv_alts[tipnode][pos]=0
            else:
                tipnode_snv_alts[tipnode]={}
                for pos in all_snvs_pos:
                    tipnode_snv_alts[tipnode][pos]=0
            local_ploidy=0
            tipnode_snv_refs[tipnode]={}
            for pos in all_snvs_pos:
                while pos>=tipnode_cnvs_pos_changes[tipnode][0][0]:
                    local_ploidy+=tipnode_cnvs_pos_changes[tipnode].pop(0)[1]
                tipnode_snv_refs[tipnode][pos]=local_ploidy-tipnode_snv_alts[tipnode][pos]
        return nodes_vars,tipnode_snv_alts,tipnode_snv_refs,tipnode_cnvs

    def tree2nhx(self,with_lens=False,attrs=None):
        '''
        Convert tree structure to string in Newick/NHX format.
        '''
        newick_str=''

        if self.left!=None:
            newick_str+='(' + self.left.tree2nhx(with_lens=with_lens,attrs=attrs)
        if self.name==None:
            newick_str+=','
        else:
            newick_str+=self.name
        if self.right!=None:
            newick_str+=self.right.tree2nhx(with_lens=with_lens,attrs=attrs) + ')'
        if self.lens!=None and with_lens:
            newick_str+= ':' + str(self.lens)
        if attrs!=None:
            newick_str+='[&&NHX'
            for attribute in attrs:
                if getattr(self, attribute):
                    if isinstance(getattr(self, attribute),list):
                        newick_str+=':{}=LIST{{{}}}'.format(attribute,'@'.join([str(x) for x in getattr(self, attribute)]))
                    elif isinstance(getattr(self, attribute),set):
                        newick_str+=':{}=SET{{{}}}'.format(attribute,'@'.join([str(x) for x in getattr(self, attribute)]))
                    elif isinstance(getattr(self, attribute),dict):
                        tmp=getattr(self, attribute)
                        newick_str+=':{}=DICT{{{}}}'.format(attribute,'@'.join(['{}>{}'.format(x,tmp[x]) for x in tmp]))
                    else:
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
    Return a list of lists. Each sublist contain two elements. The first is the position, and the second 
    is the copy number CHANGES across all the samples between that positon and the next position.
    [[pos,relative_copy_number_change],...]
    '''
    pos_changes=[[0,background],[length,-background]]
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
                mytree.lens=float(lens)
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
                raise ShouldNotBeHereError('Check your newick tree! Or maybe there are something I do not know about newick!')
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

def hap_local_leaves(positions=None,haps_cnvs=None,length=None,background=None,ploidy=None):
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
        haps_cnvs[i].sort(key=lambda cnv: cnv['start'])
        hap_cnvs_pos_changes.append(cnvs2pos_changes(cnvs=haps_cnvs[i],length=length,background=background))
        hap_local_copy.append(0)
    for pos in positions:
        all_pos_local_copy.append([pos])
        for i in range(ploidy):
            while pos>=hap_cnvs_pos_changes[i][0][0]:
                hap_local_copy[i]+=hap_cnvs_pos_changes[i].pop(0)[1]
            all_pos_local_copy[-1].append(hap_local_copy[i])
    return all_pos_local_copy

def output_tipnode_hap(tipnode_hap=None,directory=None,chroms=None,haplotype=None,parental=None):
    '''
    Output the variants of each tipnode in the order of coordinate.
    '''
    for tipnode in tipnode_hap['vars'].keys():
        with open(os.path.join(directory,'{}.genome.chain'.format(tipnode)),'a') as chain_file:
            chain_file.write('>{}_Hap{} parental:{}\n'.format(chroms,haplotype,parental))
            retrieve_tip_vars(tip_vars=tipnode_hap,tip=tipnode,out_file=chain_file,chroms=chroms)

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
                out_file.write(build_line(elements=[chroms,breakpoint,var['start'],'REF']))
                out_file.write(build_line(elements=[chroms,var['start'],var['end'],var['type'],var['mutation']]))
            elif var['start']==breakpoint:
                out_file.write(build_line(elements=[chroms,var['start'],var['end'],var['type'],var['mutation']]))
            else:
                raise ShouldNotBeHereError
            breakpoint=var['end']
        elif var['type']=='DEL': #deletion
            if var['start']>breakpoint:
                out_file.write(build_line(elements=[chroms,breakpoint,var['start'],'REF']))
                out_file.write(build_line(elements=[chroms,var['start'],var['end'],var['type'],var['copy']]))
            elif var['start']==breakpoint:
                out_file.write(build_line(elements=[chroms,var['start'],var['end'],var['type'],var['copy']]))
            else:
                raise ShouldNotBeHereError
            breakpoint=var['end']
        elif var['type']=='AMP': #amplification
            if var['start']>breakpoint:
                out_file.write(build_line(elements=[chroms,breakpoint,var['start'],'REF']))
                breakpoint=var['start']
            out_file.write(build_line(elements=[chroms,var['start'],var['end'],var['type'],'+{}'.format(var['copy'])]))
            for haplotype in var['haplotypes']:
                retrieve_tip_vars(tip_vars=haplotype,tip=tip,out_file=out_file,chroms=chroms)
        else: 
            raise ShouldNotBeHereError
    if tip_vars['end']>breakpoint:
        out_file.write(build_line(elements=[chroms,breakpoint,tip_vars['end'],'REF']))
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

class TooManyHitsOnSamePosError(Exception):
    pass

class ShouldNotBeHereError(Exception):
    pass

class TreePruneError(Exception):
    pass
