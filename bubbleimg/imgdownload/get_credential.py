# get_credential.py
# ALS 2017/05/24


import os
import getpass

def getCrendential(cred_env, cred_name = 'username'):
	cred_from_env = os.environ.get(cred_env, '')
	if cred_from_env != '':
		return cred_from_env
	else:
		print(("To avoid being asked, please set up environment variables \n     > export {cred_env} \n	> read -s {cred_env}".format(cred_env=cred_env)))

		return getpass.getpass(cred_name+': ')
