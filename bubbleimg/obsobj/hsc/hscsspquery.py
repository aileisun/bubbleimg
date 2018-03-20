""" 
do sql search on hsc data base

this is modified from hscSspQuery.py version 20160120.1 by ALS on 2017/05/11

  --------------------------------------------------------------------------------
  Original Instructions:
  https://hscdata.mtk.nao.ac.jp/hsc_ssp/dr1/common/cas_script.html
  https://hsc-gitlab.mtk.nao.ac.jp/snippets/13
  --------------------------------------------------------------------------------
  usage:

  $ echo "SELECT now();" > test.sql
  $ python hscSspQuery.py test.sql -u "your_STARS_account" > result.csv
  password? (input your STARS password)

  OR

  bash)
  ### input your STARS username and password
  $ export HSC_SSP_CAS_USERNAME
  $ read -s HSC_SSP_CAS_USERNAME
  $ export HSC_SSP_CAS_PASSWORD
  $ read -s HSC_SSP_CAS_PASSWORD

  $ python hscSspQuery.py test1.sql -u "your_STARS_account" > result1.csv
  $ python hscSspQuery.py test2.sql -u "your_STARS_account" > result2.csv
  $ python hscSspQuery.py test3.sql -u "your_STARS_account" > result3.csv
  --------------------------------------------------------------------------------

"""

import json
import argparse
import urllib.request, urllib.error, urllib.parse
import time
import sys
import csv
import getpass
import os
import os.path
import re
import ssl

version = 20160120.1


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

def hscSspQuery_retry(n_trials=20, **kwargs):
    """ 
    executing hscSspQuery and retries (n_trials times) when urllib2.URLError or urllib2.URLError happens. 

    Params
    ------
    n_trials
    **kwargs:
        arguments of hscSspQuery
    """

    for _ in range(n_trials):
        try:
            hscSspQuery(**kwargs)
            break
        except (urllib.error.HTTPError, urllib.error.URLError) as e:
            print(("[hscsspquery] retrying as error detected: "+str(e)))


def hscSspQuery(sql, filename_out='results.csv', **kwargs):

    """
    Params
    ------
    required: 
        sql (string)

    optional:
        filename_out: path of output filename
            = 'results.csv'
        release_version
            ='dr1'
        delete_job: 'delete the job you submitted after your downloading'
            =True
        out_format: choices=['csv', 'csv.gz', 'sqlite3', 'fits']
            ='csv'
        nomail: 'suppress email notice'
            =True
        username_env: 'specify the environment variable that stores STARS username'
            ='HSC_SSP_CAS_USERNAME'
        password_env: 'specify the environment variable that stores STARS password'
            ='HSC_SSP_CAS_PASSWORD'
        preview: 'quick mode (short timeout)'
            =True
        skip_syntax_check: 'skip syntax check'
            =True
        api_url
            ='https://hscdata.mtk.nao.ac.jp/datasearch/api/catalog_jobs/'

    """

    kwargs.setdefault('username_env', 'HSC_SSP_CAS_USERNAME')
    kwargs.setdefault('password_env', 'HSC_SSP_CAS_PASSWORD')
    kwargs.setdefault('release_version', 'dr1')
    kwargs.setdefault('delete_job', True)
    kwargs.setdefault('out_format', 'csv')
    kwargs.setdefault('nomail', True)
    kwargs.setdefault('preview', True)
    kwargs.setdefault('skip_syntax_check', True)

    if kwargs['release_version'][:2] == 'dr': # internal data release
        kwargs.setdefault('api_url', 'https://hscdata.mtk.nao.ac.jp/datasearch/api/catalog_jobs/')
    elif kwargs['release_version'][:3] == 'pdr': # external data release
        kwargs.setdefault('api_url', 'https://hsc-release.mtk.nao.ac.jp/datasearch/api/catalog_jobs/')
    else:
        raise ValueError("[hscsspquery] release_version not understood")


    args = Struct(**kwargs)

    credential = {'account_name': getUsername(args), 'password': getPassword(args)}
    # sql = args.__dict__['sql-file'].read()

    job = None

    try:
        if args.preview:
            with open(filename_out, 'w') as out:
                preview(credential, sql, out, args)
        else:
            job = submitJob(credential, sql, args)
            blockUntilJobFinishes(credential, job['id'], args)
            with open(filename_out, 'w') as out:
                download(credential, job['id'], out, args)
            if args.delete_job:
                deleteJob(credential, job['id'], args)
    except urllib.error.HTTPError as e:
        if e.code == 401:
            print('invalid id or password.', file=sys.stderr)
        if e.code == 406:
            print(e.read(), file=sys.stderr)
        else:
            print(e, file=sys.stderr)
    except QueryError as e:
        print(e, file=sys.stderr)
    except KeyboardInterrupt:
        if job is not None:
            jobCancel(credential, job['id'], args)
        raise KeyboardInterrupt


class QueryError(Exception):
    pass


def httpJsonPost(url, data):
    data['clientVersion'] = version
    postData = json.dumps(data)
    return httpPost(url, postData, {'Content-type': 'application/json'})


def httpPost(url, postData, headers):
    req = urllib.request.Request(url, data=postData.encode("utf-8"), headers=headers)
    skipVerifying = None
    try:
        skipVerifying = ssl.SSLContext(ssl.PROTOCOL_TLSv1)
    except AttributeError:
        pass
    if skipVerifying:
        res = urllib.request.urlopen(req, context=skipVerifying).read().decode('utf-8')
    else:
        res = urllib.request.urlopen(req).read().decode('utf-8')
    return res


def submitJob(credential, sql, args):
    url = args.api_url + 'submit'
    catalog_job = {
        'sql'                     : sql,
        'out_format'              : args.out_format,
        'include_metainfo_to_body': True,
        'release_version'         : args.release_version,
    }
    postData = {'credential': credential, 'catalog_job': catalog_job, 'nomail': args.nomail, 'skip_syntax_check': args.skip_syntax_check}
    res = httpJsonPost(url, postData)
    job = json.loads(res)
    return job


def jobStatus(credential, job_id, args):
    url = args.api_url + 'status'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    job = json.loads(res)
    return job


def jobCancel(credential, job_id, args):
    url = args.api_url + 'cancel'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


def preview(credential, sql, out, args):
    url = args.api_url + 'preview'
    catalog_job = {
        'sql'             : sql,
        'release_version' : args.release_version,
    }
    postData = {'credential': credential, 'catalog_job': catalog_job}
    res = httpJsonPost(url, postData)
    result = json.loads(res)
    writer = csv.writer(out)
    writer.writerow(result['result']['fields'])

    for row in result['result']['rows']:
        writer.writerow(row)

    if result['result']['count'] > len(result['result']['rows']):
        raise QueryError('only top %d records are displayed !' % len(result['result']['rows']))


def blockUntilJobFinishes(credential, job_id, args):
    max_interval = 5 * 60 # sec.
    interval = 1
    while True:
        time.sleep(interval)
        job = jobStatus(credential, job_id, args)
        if job['status'] == 'error':
            raise QueryError('query error: ' + job['error'])
        if job['status'] == 'done':
            break
        interval *= 2
        if interval > max_interval:
            interval = max_interval


def download(credential, job_id, out, args):
    url = args.api_url + 'download'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    bufSize = 64 * 1<<10 # 64k
    while True:
        buf = res.read(bufSize)
        out.write(buf)
        if len(buf) < bufSize:
            break


def deleteJob(credential, job_id, args):
    url = args.api_url + 'delete'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


def getPassword(args):
    password_from_envvar = os.environ.get(args.password_env, '')
    if password_from_envvar != '':
        return password_from_envvar
    else:
        return getpass.getpass('STARs password: ')


def getUsername(args):
    username_from_envvar = os.environ.get(args.username_env, '')
    if username_from_envvar != '':
        return username_from_envvar
    else:
        return getpass.getpass('STARs username: ')


if __name__ == '__main__':
    main()