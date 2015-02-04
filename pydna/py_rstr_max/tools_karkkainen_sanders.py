#! /usr/bin/python
# -*- coding: utf-8 -*-
from array import array

def radixpass(a, b, r, s, n, k) :
  c = array("i", [0]*(k+1))
  for i in xrange(n) :
    c[r[a[i]+s]]+=1

  somme = 0
  for i in xrange(k+1):
    freq, c[i] = c[i], somme
    somme += freq

  for i in xrange(n) :
    b[c[r[a[i]+s]]] = a[i]
    c[r[a[i]+s]] += 1

def simple_kark_sort(s) :
  alphabet = [None] + sorted(set(s))
  k = len(alphabet)
  t = dict((c, i) for i,c in enumerate(alphabet))
  n = len(s)
  SA = array('i', [0]*(n+3))
  s = array('i', [t[c] for c in s]+[0]*3)
  kark_sort(s, SA, n, k)
  return (s,SA)
#  return (string, SA, length, k)

def direct_kark_sort(s) :
  alphabet = [None] + sorted(set(s))
  k = len(alphabet)
  n = len(s)
  t = dict((c, i) for i,c in enumerate(alphabet))
  SA = array('i', [0]*(n+3))
  kark_sort(array('i', [t[c] for c in s]+[0]*3), SA, n, k)
  return SA[:n]

def kark_sort(s, SA, n, K):
  n0  = (n+2) / 3
  n1  = (n+1) / 3
  n2  = n / 3
  n02 = n0 + n2
      
  SA12 = array('i', [0]*(n02+3))
  SA0  = array('i', [0]*n0)

  s12 = [i for i in xrange(n+(n0-n1)) if i%3] 
  s12.extend([0]*3)
  s12 = array('i', s12)

  radixpass(s12, SA12, s, 2, n02, K)
  radixpass(SA12, s12, s, 1, n02, K)
  radixpass(s12, SA12, s, 0, n02, K)

  name = 0
  c0, c1, c2 = -1, -1, -1
  for i in xrange(n02) :
    if s[SA12[i]] != c0 or s[SA12[i]+1] != c1 or s[SA12[i]+2] != c2 :
      name += 1
      c0 = s[SA12[i]]
      c1 = s[SA12[i]+1]
      c2 = s[SA12[i]+2]
    if SA12[i] % 3 == 1 :
      s12[SA12[i]/3] = name
    else :
      s12[SA12[i]/3 + n0] = name

  if name < n02 :
    kark_sort(s12, SA12, n02, name+1)
    for i in xrange(n02) :
      s12[SA12[i]] = i+1
  else :
    for i in xrange(n02) :
      SA12[s12[i]-1] = i

  s0 = array('i',[SA12[i]*3 for i in xrange(n02) if SA12[i]<n0])
  radixpass(s0, SA0, s, 0, n0, K)
  
  p = j = k = 0
  t = n0 - n1
  while k < n :
    i = SA12[t]*3+1 if SA12[t]<n0 else (SA12[t] - n0)*3 + 2
    j = SA0[p] if p < n0 else 0

    if SA12[t] < n0 :
      test = (s12[SA12[t]+n0] <= s12[j/3]) if(s[i]==s[j]) else (s[i] < s[j])
    elif(s[i]==s[j]) :
      test = s12[SA12[t]-n0+1] <= s12[j/3 + n0] if(s[i+1]==s[j+1]) else s[i+1] < s[j+1]
    else :
      test = s[i] < s[j]

    if(test) :
      SA[k] = i
      t += 1
      if t == n02 :
        k += 1
        while p < n0 :
          SA[k] = SA0[p]
          p += 1
          k += 1
        
    else : 
      SA[k] = j
      p += 1
      if p == n0 :
        k += 1
        while t < n02 :
          SA[k] = (SA12[t] * 3) + 1 if SA12[t] < n0 else ((SA12[t] - n0) * 3) + 2
          t += 1
          k += 1
    k += 1

def LCP(s, suffix_array):
  n = len(s)
  init = [0]*n
  rank = array('i', init)
  LCP = array('i', init)
  for i in xrange(n):
    rank[suffix_array[i]] = i
  l = 0
  for j in xrange(n):
    l = max(0, l-1)
    i = rank[j]
    j2 = suffix_array[i-1]
    if i:
      while l + j < n and l + j2 < n and s[j+l] == s[j2+l]:
        l += 1
      LCP[i-1] = l
    else:
      l = 0
  return LCP

