#!/usr/bin/python
#-*- coding: utf-8 -*-

#   This software was developed by Zanoni Dias, Ulisses Dias and João
#   C. Setubal
#
#   It should not be redistributed or used for any commercial purpose
#   without written permission from authors
#
#   release date: nov 15, 2011
#
# This software is experimental in nature and is
# supplied "AS IS", without obligation by the authors to provide
# accompanying services or support.  The entire risk as to the quality
# and performance of the Software is with you. The authors
# EXPRESSLY DISCLAIM ANY AND ALL WARRANTIES REGARDING THE SOFTWARE,
# WHETHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES
# PERTAINING TO MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
#

import sys
import re

# Returns the second element of an array
def second_element(array) :
    return array[1]

# Returns the least element between the first and the second element
# of an array
def min_first_second_element(array) :
    return min(array[0], array[1])

# Parse the coordinate file generated by nucmer or promer
def parse_coords(coords_filename) :
    info_array  = []
    coords_file = open(coords_filename)
    coords_file.readline()
    coords_file.readline()
    coords_file.readline()
    coords_file.readline()
    coords_file.readline()
    for line in coords_file :
        match =  re.match("\s+(\S+)\s+(\S+)\s+\|\s+(\S+)\s+(\S+)\s+\|(.*)\s(\S+)", line)
        if match :
            s1     = int(match.group(1))
            e1     = int(match.group(2))
            s2     = int(match.group(3))
            e2     = int(match.group(4))
            extra  = match.group(5)
            contig = match.group(6)
            info_array.append([s1,e1,s2,e2,contig])
    coords_file.close()
    info_array.sort(key=min_first_second_element)
    return info_array

# Receives an array of coordinates and generates the contig
# representation as a set of integers. The assumption we make
# here is that the reference is the identity permutation.
def coords_to_permutation(coords) :
    contigs         = {}
    best_scored_mum = {}
    count           = 1

    def get_orientation(lcount,ls1,le1,ls2,le2) :
        if ls1 < le1 :
            if ls2 < le2 :
                return [ lcount, ls2, le2]
            else :
                return [-lcount, le2, ls2]
        else :
            if ls2 < le2 :
                return [-lcount, ls2, le2]
            else :
                return [ lcount, le2, ls2]

    for el in coords :
        s1     = el[0]
        e1     = el[1]
        s2     = el[2]
        e2     = el[3]
        contig = el[4]
        aux = get_orientation(count, s1, e1, s2, e2)
        if not contig in contigs :
            contigs[contig] = []
        contigs[contig].append([aux[0], aux[1], aux[2], contig])
        count = count + 1

    contigs_array = []
    for contig in contigs.keys() :
        contigs[contig].sort(key=second_element)
        contigs_array.append([contig])
        for el in contigs[contig] :
            contigs_array[-1].append(el[0])

    return contigs_array


##################################################################
class InContigElement:
  def __init__(self,contig,position,sign):
    self.contig = contig
    self.position = position
    self.sign = sign

##################################################################
# Return the sign of an element
def sign(n):

  return n / abs(n)

##################################################################
# Return an array of size n filled with zeros
def zeros(n):

  z = []
  for i in range(n):
    z.append(0)

  return z

##################################################################
# Create the data structure that allows the search for the contig
# where one element is placed in.
def makeInContig(contigs,c):

  inContig = []
  for i in range(n+1):
    inContig.append(None)

  for i in range(c):
    for j in range(len(contigs[i])):
      element = abs(contigs[i][j])
      inContig[element] = InContigElement(i,j,sign(contigs[i][j]))

  return inContig

##################################################################
# Generate the scaffold from a set of contigs. This procedure is the
# main core of our method.
def makeScaffold2(contigs,inContig,c,n):
  scaffold = []
  lookingfor = 1
  used = zeros(c)
  next = 0
  s = 0
  fail = 0
  ok = 0
  searched = []
  restart    = True
  breakpoint = True

  # While there is contig to place in the righ position
  while s < c:
    reset = 0

    # If we are looking for a valid element
    if lookingfor != n+1:

      element = inContig[abs(lookingfor)]

      if sign(lookingfor) == element.sign:
        if element.position != 0:
          lookingfor = contigs[element.contig][element.position-1] + 1
          breakpoint = True
        else:
          ok = 1
      else:
        if element.position != len(contigs[element.contig])-1:
          lookingfor = - (contigs[element.contig][element.position+1] - 1)
          breakpoint = True
        else:
          ok = 1

    # If we are looking for a new element
    if searched.count(lookingfor) == 0:
      searched.append(lookingfor)
    else:
      ok = 1

    if ok:
      if not(used[element.contig]):
        scaffold.append([sign(lookingfor) * (element.contig+1) * element.sign,restart,breakpoint])
        restart    = False
        breakpoint = False
        s += 1
        searched = []
        used[element.contig] = 1
        while next < c and used[next]:
          next += 1

        if sign(lookingfor) == element.sign:
          lookingfor = contigs[element.contig][-1] + 1
        else:
          lookingfor = -contigs[element.contig][0] + 1

        ok = 0
      else:
        reset = 1

    # Check if it is necessary to start a new search using a
    # non-placed contig
    if (s < c) and (lookingfor == n+1 or reset):
      lookingfor = contigs[next][0]
      restart    = True
      breakpoint = True
      searched = []
      fail += 1

  return [scaffold, fail]

#######################################################
####################### Main ##########################
#######################################################

try :
    # Parse the input
    coords = parse_coords(sys.argv[1])

    # Generate the contig array
    contig_array = coords_to_permutation(coords)
    contigs = []
    for contig in contig_array :
        contigs.append(contig[1:])

    c = len(contigs)
    n = len(coords)

    # Generate the inContig structure. This structure make it easier
    # to search for the contig where a given element is placed in
    inContig = makeInContig(contigs,c)

    # Generate the scaffold
    scaffold = makeScaffold2(contigs, inContig, c, n)

    # Print the scaffold

    scaffold_count   = 1
    print("> Scaffold_%s" % scaffold_count)
    for el in scaffold[0] :
        if el[1] :
            if scaffold_count > 1 :
                print("")
                print("> Scaffold_%s" % scaffold_count)
            scaffold_count = scaffold_count + 1

        if el[0] < 0 :
            print("%s 1" % contig_array[int(abs(el[0])-1)][0])
        else :
            print("%s 0" % contig_array[int(el[0]-1)][0])
except :
    print("Usage: %s <coord file>" % sys.argv[0])
    print("   <coord file> is the coordinate file of show-coord (MUMmer)")
