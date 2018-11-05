#!/usr/bin/env python

#########################################################################
# Author: Hechuan Yang
# Created Time: 2016-12-21 18:00:34
# Description: This script is used to process the output of tumopp
#########################################################################

import os
import sys
import re
import argparse
import copy
import fileinput

#handle the error below
#python | head == IOError: [Errno 32] Broken pipe 
from signal import signal, SIGPIPE, SIG_DFL 
signal(SIGPIPE,SIG_DFL) 

class tree:
    def __init__(self,name=None,lens=None,left=None,right=None,top=None,
        nodeid=None,x=None,y=None,z=None,birth=None,death=None,omega=None):
        self.name=name
        self.lens=lens
        self.left=left
        self.right=right
        self.top=top
        self.nodeid=nodeid
        self.x=x
        self.y=y
        self.z=z
        self.birth=birth
        self.death=death
        self.omega=omega

    def add_node(self,node):
        '''
        This method is used when constructing a tree.
        After adding node to the growing tree, return the new node added.
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

    def ms_labeling(self):
        '''
        Label all the leaves with 1,2,3,...,
        and then label the inner nodes with the smallest one of the leaves' labels it has.
        '''
        all_leaves=self.leaves_naming()
        try:
            all_leaves=[str(x) for x in sorted([int(x) for x in all_leaves])]
        except ValueError:
            all_leaves.sort()
        label_map={}
        for i in range(len(all_leaves)):
            label_map[all_leaves[i]]=i+1
        self.translate_name(label_map=label_map)
        return label_map

    def translate_name(self,label_map=None):
        '''
        Translate leaf name to ms_label.
        '''
        if not hasattr(self,'ms_label') or self.ms_label == None:
            if self.left==None and self.right==None:
                self.ms_label=label_map[self.name]
            else:
                self.left.translate_name(label_map=label_map)
                self.right.translate_name(label_map=label_map)
                self.ms_label=str(sorted([int(self.left.ms_label),int(self.right.ms_label)])[0])
        return self.ms_label

    def ms_ejoin_collecting(self,now=None):
        '''
        Generate the ms populaiton split string of the the tree.
        '''
        if not hasattr(self,'ms_ejoin') or self.ms_ejoin == None:
            if self.left==None and self.right==None:
                self.ms_ejoin=''
            else:
                self.ms_ejoin=self.left.ms_ejoin_collecting(now=now)+self.right.ms_ejoin_collecting(now=now)+\
                    ' -ej {} '.format(str(now-self.death))+\
                    ' '.join([str(x) for x in reversed(sorted([int(self.right.ms_label),int(self.left.ms_label)]))])
        return self.ms_ejoin
    
    def nhx_tree(self,lens=False,attr=None):
        nhx_str=self.nhx(lens=lens,attr=attr)
        if nhx_str:
            nhx_str+=';'
        return nhx_str

    def nhx(self,lens=False,attr=None):
        nhx_str=''
        if self is None: 
            return nhx_str
        if self.left != None:
            nhx_str+='('+self.left.nhx(lens=lens,attr=attr)
        if self.name==None:
            nhx_str+=','
        else:
            nhx_str+=self.name

        if self.right != None:
            nhx_str+=self.right.nhx(lens=lens,attr=attr)+')'
        if self.lens!=None and lens:
            nhx_str+=':{}'.format(self.lens)
        if attr!=None :
            nhx_str+='[&&NHX'
            for attribute in attr:
                if getattr(self, attribute,None)!=None:
                    nhx_str+=':{}={}'.format(attribute,getattr(self, attribute))
            nhx_str+=']'
        return nhx_str

    def subtree(self,selected=None):
        '''
        Return the subtree of the tipnodes in selcted.
        '''
        assert isinstance(selected,set), \
            'The selected should be a set containing the names of tipnodes!'
        subset=selected.intersection(self.leaves_naming())
        if subset:
            subtree1=None
            subtree2=None
            newnode=tree(
                name=self.name,
                lens=self.lens,
                left=None,
                right=None,
                top=self.top,
                nodeid=self.nodeid,
                x=self.x,
                y=self.y,
                z=self.z,
                birth=self.birth,
                death=self.death,
                omega=self.omega)
            if hasattr(self,'sector'):
                newnode.sector=self.sector

            if self.left!=None:
                subtree1=self.left.subtree(selected=subset)
            if self.right!=None:
                subtree2=self.right.subtree(selected=subset)
            if subtree1 and subtree2:
                newnode.left=subtree1
                subtree1.top=newnode
                newnode.right=subtree2
                subtree2.top=newnode
            elif subtree1:
                subtree1.top=newnode.top
                if subtree1.lens!=None:
                    subtree1.lens+=newnode.lens
                if newnode.top!=None:
                    newnode.top.death=subtree1.birth
                    if newnode.top.left is newnode:
                        newnode.top.left=subtree1
                    if newnode.top.right is newnode:
                        newnode.top.right=subtree1
                newnode=subtree1
            elif subtree2:
                subtree2.top=newnode.top
                if subtree2.lens!=None:
                    subtree2.lens+=newnode.lens
                if newnode.top!=None:
                    newnode.top.death=subtree2.birth
                    if newnode.top.left is newnode:
                        newnode.top.left=subtree2
                    elif newnode.top.right is newnode:
                        newnode.top.right=subtree2
                newnode=subtree2
            return newnode
        else:
            return None

    def tree2ms(self,now=None,psam=None,howmany=None):
        '''
        I will generate the ms code of the history of the tree.
        '''
        nsam=psam*self.leaves_counting()
        mscode='ms {} {} -T -I {}'.format(nsam,howmany,self.leaves_counting())
        for i in range(self.leaves_counting()):
            mscode+=' {}'.format(psam)
        label_map=self.ms_labeling()
        mscode+=self.ms_ejoin_collecting(now=now)
        return label_map,mscode
            
def newick2tree(newick=None):
    leaf_name_re=re.compile('\w+:')
    lens_re=re.compile(':[0-9.]+')
    while newick != ';':
        if newick.startswith('('):
            if 'mytree' in vars():
                mytree=mytree.add_node(tree())
            else:
                mytree=tree()
            newick=newick[1:]
        elif leaf_name_re.match(newick):
            m=leaf_name_re.match(newick)
            index=m.span()
            leaf_name=newick[:index[1]-1]
            newick=newick[index[1]-1:]
            mytree=mytree.add_node(tree(name=leaf_name))
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
            print('Should not be here!2')
            print('Check your newick tree! Or maybe there are something I do not know about newick!')
    return mytree

def tumopp2tree(tumopp=None,sectors=None,justpick=None):
    '''
    Build tree from tumopp output, and output three variables:
    root:
    nodes:
    picked:
    limits:
    '''
    largest = float('Inf')
    smallest = -float('Inf')
    root=None
    nodes={}
    picked={}
    for sector in sectors:
        picked[sector]=[]
    limits={'x':[largest,smallest],
            'y':[largest,smallest],
            'z':[largest,smallest]}
              
    input=tumopp
    header=next(input)
    header=header.rstrip()
    header=header.split()
    for line in input:
        line=line.rstrip()
        columns=line.split()
        info={}
        for tag,value in zip(header,columns):
            if tag in ('id','ancestor'):
                info[tag]=value
            else:
                info[tag]=float(value)

        if info['death']==0:
            if limits['x'][0]>info['x']:
                limits['x'][0]=info['x']
            if limits['x'][1]<info['x']:
                limits['x'][1]=info['x']
            if limits['y'][0]>info['y']:
                limits['y'][0]=info['y']
            if limits['y'][1]<info['y']:
                limits['y'][1]=info['y']
            if limits['z'][0]>info['z']:
                limits['z'][0]=info['z']
            if limits['z'][1]<info['z']:
                limits['z'][1]=info['z']
                
            for sector in sectors:
                xrange,yrange,zrange=sectors[sector]
                if xrange[0]<info['x']<xrange[1] and \
                    yrange[0]<info['y']<yrange[1] and \
                    zrange[0]<info['z']<zrange[1]:
                    picked[sector].append(info['id'])
                    info['sector']=sector
                    break
        if justpick:
            continue

        node=tree(nodeid=info['id'],
            x=info['x'],
            y=info['y'],
            z=info['z'],
            birth=info['birth'],
            death=info['death'],
            omega=info['omega'])
        if node.death!=0:
            node.lens=node.death-node.birth
        else:
            node.name=info['id']
            if 'sector' in info:
                node.sector=info['sector']

        if node.nodeid=='1':
            root=node

        if node.nodeid in nodes:
            if nodes[node.nodeid].left!=None:
                node.left=nodes[node.nodeid].left
                nodes[node.nodeid].left.top=node
            if nodes[node.nodeid].right!=None:
                node.right=nodes[node.nodeid].right
                nodes[node.nodeid].right.top=node
        nodes[node.nodeid]=node

        if info['ancestor'] in nodes:
            if nodes[info['ancestor']].left==None:
                nodes[info['ancestor']].left=node
                node.top=nodes[info['ancestor']]
            elif nodes[info['ancestor']].right==None:
                nodes[info['ancestor']].right=node
                node.top=nodes[info['ancestor']]
            if nodes[info['ancestor']].left and nodes[info['ancestor']].right:
                if int(nodes[info['ancestor']].left.nodeid)>int(nodes[info['ancestor']].right.nodeid):
                    nodes[info['ancestor']].left,nodes[info['ancestor']].right=\
                    nodes[info['ancestor']].right,nodes[info['ancestor']].left
        elif info['ancestor']!='0':
            nodes[info['ancestor']]=tree(nodeid=info['ancestor'],left=node)
            node.top=nodes[info['ancestor']]
    return root,nodes,picked,limits

def assign_tipnode_lens(nodes=None,shortest=None):
    '''
    In tumopp output, the alive cells' death is 0.
    In this case, the length of all tipnode is None.
    I will assign a length to these tipnodes. You can specify the shortest length,
    which will be assigned to the youngest cell. Otherwise, I will pick the shortest
    leangth of all non-tipnode as the shortest length.
    I will also return the time of 'now'
    '''
    if shortest==None:
        all_lens=[node.lens for node in nodes.values() if node.lens]
        shortest=all_lens[0]
        for lens in all_lens:
            if lens<shortest:
                shortest=lens
    latest_birth=0
    for node in nodes.values():
        if node.birth>latest_birth:
            latest_birth=node.birth
    now=latest_birth+shortest
    for node in nodes.values():
        if node.death==0:
            node.death=now
            node.lens=node.death-node.birth
    return now

def check_range(s=None):
    try:
        start,end=s.split(',')
        start=float(start)
        end=float(end)
    except ValueError:
        argparse.ArgumentTypeError('{} is not a valid value for range!'.format(s))
    return start,end

def parse_sector(sectorf=None):
    '''
    Parse the range file into a dictionary with the structure:
    {sector1:[[xstart,xend],[ystart,yend],[zstart,zend]],
     sector2:[[xstart,xend],[ystart,yend],[zstart,zend]],
     ...
    }
    If the value of one coordinate range of a sector is ',', 
    will just put None there instead.
    '''
    sectors={}
    with open(sectorf) as input:
        header=next(input)
        header=header.rstrip()
        header=header.split()
        assert set(header)==set(['sector','x','y','z']), \
            "The range file should start with a header containing 'sector','x','y','z'!"
        for line in input:
            line=line.rstrip()
            columns=line.split()
            info={}
            for tag,value in zip(header,columns):
                info[tag]=value
            sectors[info['sector']]=[]
            for coordinate in 'x','y','z':
                if info[coordinate]==',':
                    sectors[info['sector']].append(None)
                else:
                    start,end=info[coordinate].split(',')
                    start=float(start)
                    end=float(end)
                    sectors[info['sector']].append([start,end])
    return sectors

def parse_group(groupf=None):
    '''
    Read the group settings from group file, which contains two columns with a header line:
    sector: sector name
    id: the nodeids in the sector separated by commas
    '''
    sectors={}
    with open(groupf) as input:
        header=next(input)
        header=header.rstrip()
        header=header.split()
        assert set(header)==set(['sector','id']), \
            "The range file should start with a header containing 'sector','id'!"
        for line in input:
            line=line.rstrip()
            columns=line.split()
            sectors[columns[0]]=columns[1].split(',')
    return sectors

def parse_range(sectors=None):
    '''
    Convert the None value in range coordinate into [-Inf,Inf].
    '''
    largest = float('Inf')
    smallest = -float('Inf')
    for sector in sectors:
        xrange,yrange,zrange=sectors[sector]
        if xrange==None:
            xrange=[smallest,largest]
        if yrange==None:
            yrange=[smallest,largest]
        if zrange==None:
            zrange=[smallest,largest]
        sectors[sector]=[xrange,yrange,zrange]
    return sectors

#In this program, we define the population in three leveles:
#sector, deme and cell
#1. a cell is a cell. 
#2. we treat a record of tumopp as a deme, each deme can contains 
#   a number of cells
#3. a sector is a number of demes locate in the same place
if __name__ == '__main__':
    parser=argparse.ArgumentParser(
        description='Build ms code for the sectors spatially picked from tumopp simulation')
    parser.add_argument('tumopp',metavar='INPUT', nargs=1,help="input file generated by tumopp, if '-', stdin is used")
    parser.add_argument('-x','--xrange',type=check_range,help='the range of x to pick')
    parser.add_argument('-y','--yrange',type=check_range,help='the range of y to pick')
    parser.add_argument('-z','--zrange',type=check_range,help='the range of z to pick')
    default=None
    parser.add_argument('-s','--sector',type=str,default=default,
        help='the input file which contains the spatial ranges of sectors [{}]'.format(default))
    default=None
    parser.add_argument('-g','--group',type=str,default=default,
        help='the input file which contains the IDs of nodes in sectors [{}]'.format(default))
    default=None
    parser.add_argument('-T','--tree',type=str,default=default,
        help='the output file to save the original tree simulated by tumopp [{}]'.format(default))
    default=None
    parser.add_argument('-t','--subtree',type=str,default=default,
        help='the output file to save the subtree picked by coordinate [{}]'.format(default))
    default='sampled.txt'
    parser.add_argument('-S','--sample',type=str,default=default,
        help='the output file to save the samples picked by coordinate [{}]'.format(default))
    default=1
    parser.add_argument('-p','--psam',type=int,default=default,
        help='the number of samples from each deme [{}]'.format(default))
    default=1
    parser.add_argument('-n','--howmany',type=int,default=default,
        help='the number of independent samples to generate [{}]'.format(default))
    default=None
    parser.add_argument('-m','--map',type=str,default=default,
        help='the output file to save the map between tumopp deme id and ms population id [{}]'.format(default))
    default=None
    parser.add_argument('--ms',type=str,default=default,
        help='the output file to save the ms code [{}]'.format(default))
    parser.add_argument('--justpick',action="store_true",
        help='just output the slected node records'.format(default))
    parser.add_argument('--summary',action="store_true",
        help='just print out the size summary of the input tumor'.format(default))
    args=parser.parse_args()

    sectors={}
    if args.sector:
        sectors=parse_sector(args.sector)
    elif args.xrange or args.yrange or args.zrange:
        sectors={'sector0':[args.xrange,args.yrange,args.zrange]}
    #print(sectors)
    sectors=parse_range(sectors)
    #print(sectors)

    if args.tumopp:
        inputf=args.tumopp
    elif '-' in sys.argv:
        inputf=args.tumopp
    else:
        parser.print_help()
        exit(0)
    with fileinput.input(files=inputf) as tumopp:
        if args.group:
#I will overide other sector settings, if user specify the sectors explicitly by --group
            sectors={}
            root,nodes,picked,limits=tumopp2tree(tumopp=tumopp,sectors=sectors,justpick=args.justpick)
            picked=parse_group(args.group)
            for sector in picked:
                for node in picked[sector]:
                    assert node in nodes and nodes[node].death==0, \
                        '{} in sector {} is not a tipnode on the Tumopp tree!'.format(node,sector)
                    nodes[node].sector=sector
        else:
            root,nodes,picked,limits=tumopp2tree(tumopp=tumopp,sectors=sectors,justpick=args.justpick)

    #print(picked)
    if args.summary:
        print('Here is the size summary of the input tumor:')
        print('xlimits:\t'+'\t'.join([str(i) for i in limits['x']]))
        print('ylimits:\t'+'\t'.join([str(i) for i in limits['y']]))
        print('zlimits:\t'+'\t'.join([str(i) for i in limits['z']]))
    if args.sample:
        with open(args.sample,'w') as output:
            output.write('#sector\tcount\tid\n')
            for sector in sorted(picked.keys()):
                output.write('{}\t{}\t{}\n'.format(sector,len(picked[sector]),','.join(picked[sector])))
    if args.justpick:
        exit(0)

    now=assign_tipnode_lens(nodes=nodes)
    if args.tree:
        with open(args.tree,'w') as output:
            output.write(root.nhx_tree(attr=['nodeid','birth','death'])+'\n')

    sampled=[]
    for ids in picked.values():
        sampled.extend(ids)
    #print(root.nhx_tree(attr=['nodeid','sector']))
    if sampled:
        subtree=root.subtree(selected=set(sampled))
    else:
        subtree=root
    if args.subtree:
        with open(args.subtree,'w') as output:
            output.write(subtree.nhx_tree(lens=True,attr=['nodeid','sector','lens'])+'\n')

    if args.map or args.ms:
        label_map,ms=subtree.tree2ms(now=now,psam=args.psam,howmany=args.howmany)
        if args.map:
            with open(args.map,'w') as output:
                output.write('#tumopp\tms\n')
                for deme in sorted(label_map.keys(),key=lambda deme: label_map[deme]):
                    output.write('{}\t{}\n'.format(deme,label_map[deme]))
        if args.ms:
            with open(args.ms,'w') as output:
                output.write(ms+'\n')
