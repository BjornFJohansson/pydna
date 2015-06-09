#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Provides a class for downloading sequences from genbank.
'''
import cPickle
import shelve

import re
import os
import urllib2
import warnings
import sys

from urlparse               import urlparse
from urlparse               import urlunparse
from Bio                    import Entrez
from Bio.SeqUtils.CheckSum  import seguid

from pydna.dsdna import read


from pydna.pretty import pretty_str
from pydna.dsdna import Dseqrecord

def _get_proxy_from_global_settings():
    """Get proxy settings from linux/gnome"""
    if sys.platform.startswith('linux'):
        try:
            from gi.repository import Gio
        except ImportError:
            return ''
        mode = Gio.Settings.new('org.gnome.system.proxy').get_string('mode')
        if mode == 'none' or mode == 'auto':
            return None
        http_settings = Gio.Settings.new('org.gnome.system.proxy.http')
        host = http_settings.get_string('host')
        port = http_settings.get_int('port')
        if http_settings.get_boolean('use-authentication'):
            username = http_settings.get_string('authentication_user')
            password = http_settings.get_string('authentication_password')
        else:
            username = password = None
            return 'http://{}:{}'.format(host, port)
    return None

class GenbankRecord(Dseqrecord):

    def __init__(self, acc=None, *args, **kwargs):
        super(Dseqrecord, self).__init__(self, *args, **kwargs)
        self.acc = acc

    def url(self):
        return pretty_str("http://www.ncbi.nlm.nih.gov/nucleotide/"+self.acc)



class Genbank():
    '''Class to facilitate download from genbank.

    Parameters
    ----------
    users_email : string
        Has to be a valid email address. You should always tell
        Genbanks who you are, so that they can contact you.
    proxy : string, optional
        String containing a proxy url:
        "proxy = "http://umiho.proxy.com:3128"
    tool : string, optional
        Default is "pydna". This is to tell Genbank which tool you are
        using.

    Examples
    --------

    >>> import pydna
    >>> gb=pydna.Genbank("me@mail.se", proxy = "http://proxy.com:3128")                  #doctest: +SKIP
    >>> rec = gb.nucleotide("L09137") # <- pUC19 from genbank                            #doctest: +SKIP
    >>> print len(rec)                                                                   #doctest: +SKIP
    2686                                                                                 #doctest: +SKIP
    '''

    def __init__(self, users_email, proxy = None, tool="pydna"):

        if not re.match("[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}",users_email,re.IGNORECASE):
            raise(ValueError("Not a valid email address!"))

        self.email=users_email #Always tell NCBI who you are

        #print "#####", proxy

        if proxy:
            parsed = urlparse(proxy)
            scheme = parsed.scheme
            hostname = parsed.hostname
            test = urlunparse((scheme, hostname,'','','','',))
            try:
                response=urllib2.urlopen(test, timeout=1)
            except urllib2.URLError as err:
                warnings.warn("could not contact proxy server")
            #print { scheme : parsed.geturl() }
            self.proxy = urllib2.ProxyHandler({ scheme : parsed.geturl() })
        else:
            pass
            #proxy_handler = urllib2.ProxyHandler({})
            #opener = urllib2.build_opener(proxy_handler)
            #urllib2.install_opener(opener)
            #os.environ['http_proxy']=''
            #self.proxy = urllib2.ProxyHandler()
        #self.opener = urllib2.urlopen #build_opener(self.proxy)
        #urllib2.install_opener(self.opener)

    def __repr__(self):
        return "Genbank({})".format(self.email)

    def nucleotide(self, item, start=None, stop=None, strand="watson" ):
        '''Download a genbank nuclotide record.

        Item is a string containing one genbank acession number [#]_
        for a nucleotide file. Start and stop are intervals to be
        downloaded. This is useful as some genbank records are large.
        If strand is "c", "C", "crick", "Crick", "antisense","Antisense",
        "2" or 2, the watson(antisense) strand is returned, otherwise
        the sense strand is returned.

        Alternatively, item can be a string containing an url that returns a
        sequence in genbank or FASTA format.

        The genbank E-utilities can provide such urls [#]_.

        Result is returned as a Dseqrecord object.

        Genbank nucleotide accession numbers have this format:

        | A12345   = 1 letter  + 5 numerals
        | AB123456 = 2 letters + 6 numerals


        References
        ----------

        .. [#]   http://www.dsimb.inserm.fr/~fuchs/M2BI/AnalSeq/Annexes/Sequences/Accession_Numbers.htm
        .. [#]   http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
        '''
        cached  = False
        refresh = False

        if os.environ["pydna_cache"] in ("compare", "cached"):
            cache = shelve.open(os.path.join(os.environ["pydna_data_dir"],"genbank"), protocol=cPickle.HIGHEST_PROTOCOL, writeback=False)
            key = item+str(start)+str(stop)+str(strand)
            try:
                cached = cache[key]
            except:
                if os.environ["pydna_cache"] == "compare":
                    raise Exception("no result for this key!")
                else:
                    refresh = True

        if refresh or os.environ["pydna_cache"] in ("compare", "refresh", "nocache"):

            if item.lower().startswith(("ftp","http", "https")):
                try:
                    url = urlparse(item)
                except (AttributeError, TypeError):
                    url = urlparse("")

                if url.scheme:
                    return read( urllib2.urlopen(item).read() )

            matches =((1, re.search("(REGION:\s(?P<start>\d+)\.\.(?P<stop>\d+))", item)),
                      (2, re.search("(REGION: complement\((?P<start>\d+)\.\.(?P<stop>\d+)\))",item)),
                      (1, re.search(":(?P<start>\d+)-(?P<stop>\d+)",item)),
                      (2, re.search(":c(?P<start>\d+)-(?P<stop>\d+)",item)),
                      (0, None))

            # BK006936.2 REGION: complement(613900..615202)
            # NM_005546 REGION: 1..100
            # NM_005546 REGION: complement(1..100)
            # 21614549:1-100
            # 21614549:c100-1

            # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=BK006936.2&strand=2&seq_start=613900&seq_stop=615202&rettype=fasta&retmode=text
            # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=21614549&strand=1&seq_start=1&seq_stop=100&rettype=fasta&retmode=text
            # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=21614549&strand=2&seq_start=1&seq_stop=100&rettype=fasta&retmode=text
            # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=21614549&strand=1&seq_start=1&seq_stop=100&rettype=gb&retmode=text
            # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=21614549&strand=2&seq_start=1&seq_stop=100&rettype=gb&retmode=text


            for strand, match in matches:
                if match:
                    start = match.group("start")
                    stop  = match.group("stop")
                    item = item[:match.start()]
                    break

            #        if not strand:
            #            raise Exception("\nACESSION string is malformed!\n"
            #                              "NM_005546 REGION: 1..100\n"
            #                              "NM_005546 REGION: complement(1..100)\n"
            #                              "21614549:1-100\n"
            #                              "21614549:c100-1\n")

            if str(strand).lower() in ("c","crick", "antisense", "2"):
                strand = 2
            else:
                strand = 1

            Entrez.email = self.email

            #        print item
            #        print start
            #        print stop
            #        print strand

            result = read(Entrez.efetch(db        ="nucleotide",
                                        id        = item,
                                        rettype   = "gb",
                                        seq_start = start,
                                        seq_stop  = stop,
                                        strand    = strand,
                                        retmode   = "text").read())

        if os.environ["pydna_cache"] == "compare":
            if result!=cached:
                module_logger.warning('download error')

        if refresh or os.environ["pydna_cache"] == "refresh":
            cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "genbank"), protocol=cPickle.HIGHEST_PROTOCOL, writeback=False)
            cache[key] = result

        elif cached and os.environ["pydna_cache"] not in ("nocache", "refresh"):
            result = cached
            cache.close()

        return result

class web():

    def __init__(self, proxy = None):
        if proxy:
            parsed = urlparse(proxy)
            scheme = parsed.scheme
            hostname = parsed.hostname
            test = urlunparse((scheme, hostname,'','','','',))
            try:
                response=urllib2.urlopen(test, timeout=1)
            except urllib2.URLError as err:
                warnings.warn("could not contact proxy server")
            self.proxy = urllib2.ProxyHandler({ scheme : parsed.geturl() })
        else:
            pass
            #proxy_handler = urllib2.ProxyHandler({})
            #opener = urllib2.build_opener(proxy_handler)
            #urllib2.install_opener(opener)
            #os.environ['http_proxy']=''
            #self.proxy = urllib2.ProxyHandler()
        #self.opener = urllib2.urlopen #build_opener(self.proxy)
        #urllib2.install_opener(self.opener)


    def download(self, url):
        cached  = False
        refresh = False
        cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "web"), protocol=cPickle.HIGHEST_PROTOCOL, writeback=False)
        key = str(url)

        if os.environ["pydna_cache"] in ("compare", "cached"):
            try:
                cached = cache[key]
            except KeyError:
                if os.environ["pydna_cache"] == "compare":
                    raise Exception("no result for this key!")
                else:
                    refresh = True

        if refresh or os.environ["pydna_cache"] in ("compare", "refresh", "nocache"):
            response = urllib2.urlopen(url)
            result = response.read()

        if os.environ["pydna_cache"] == "compare":
            if result!=cached:
                module_logger.warning('download error')

        if refresh or os.environ["pydna_cache"] == "refresh":
            cache = shelve.open(os.path.join(os.environ["pydna_data_dir"],"genbank"), protocol=cPickle.HIGHEST_PROTOCOL, writeback=False)
            cache[key] = result

        elif cached and os.environ["pydna_cache"] not in ("nocache", "refresh"):
            result = cached

        cache.close()

        return result

if __name__=="__main__":
    import doctest
    doctest.testmod()
