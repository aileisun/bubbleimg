import abc
import paramiko
import requests
import os

from get_credential import getCrendential

class multiButler(object):

	def __init__(self, environment='iaa', **kwargs):
		"""
		decide which subclass to use
		"""
		if environment == 'iaa':
			self.butler = iaaButler(**kwargs)
		elif environment == 'online':
			self.butler = onlineButler(**kwargs)
		elif environment == 'sumire':
			self.butler = sumireButler(**kwargs)
		else: 
			raise InputError('arg environment not recognized')


class parentButler(object):
	"""
	to access hsc files in

	 different environments (online, iaa, sumire), etc. It does not nessisarily use the lsst butler. 
	"""

	def __init__(self, release_version='dr1', semester='s16a', rerun='s16a_wide'):
		self.release_version = release_version
		self.semester = semester
		self.rerun = rerun
		self.rerun_path = self._get_rerun_path()


	@abc.abstractmethod
	def download_file(self, localpath, tract, patch_s, filter, coadd='deepCoadd', filetype='calexp'): 
		"""
		download specified file to localpath

		Params
		------
		localpath (string): path/file name to save to 
		tract (int): tract number
		patch_s (str): e.g., '7,3' patch_s number
		filter (str): e.g., 'HSC-G'
		coadd='deepCoadd'
		filetype='calexp'

		Return
		------
		status: True if successful, false if not.
		"""
		raise NotImplementedError("Subclass must implement abstract method")


	@abc.abstractmethod
	def _get_root_path(self):
		raise NotImplementedError("Subclass must implement abstract method")


	def _get_rerun_path(self):
		"""
		e.g., "dr1/s16a/data/s16a_wide/"
		"""
		return '{dr}/{semester}/data/{rerun}/'.format(dr=self.release_version, semester=self.semester, rerun=self.rerun)


	def _get_tail_path(self, tract, patch_s, filter, coadd='deepCoadd', filetype='calexp'):
		"""
		e.g., "deepCoadd/HSC-R/9564/7,3/calexp-HSC-R-9564-7,3.fits"
		"""
		return '{coadd}/{filter}/{tract}/{patch_s}/{filetype}-{filter}-{tract}-{patch_s}.fits'.format(coadd=coadd, filetype=filetype, tract=str(tract), patch_s=str(patch_s), filter=filter)



class iaaButler(parentButler):

	"""	

	Instruction
	-----------
	One needs to set sumire username and password as environment parameters

	export IAA_SUMIRE_CLUSTER_USERNAME
	read -s IAA_SUMIRE_CLUSTER_USERNAME
	export IAA_SUMIRE_CLUSTER_PASSWORD
	read -s IAA_SUMIRE_CLUSTER_PASSWORD

	"""
	def __init__(self, **kwargs):

		super(self.__class__, self).__init__(**kwargs)

		self.__username = getCrendential("IAA_SUMIRE_CLUSTER_USERNAME", cred_name = 'iaa username')
		self.__password = getCrendential("IAA_SUMIRE_CLUSTER_PASSWORD", cred_name = 'iaa password')

		self.root_path = self._get_root_path()


	def __enter__(self):
		""" Connect """
		self.ssh_client = paramiko.SSHClient()
		self.ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		self.ssh_client.connect('sumire.asiaa.sinica.edu.tw', username=self.__username, password=self.__password, port=20001)


	def __exit__(self, exc_type, exc_value, traceback):
		""" Close """
		self.ssh_client.close()


	def download(self, remotepath, localpath):
		""" force download """

		if os.path.isfile(localpath):
			os.remove(localpath)

		try:
			with self.ssh_client.open_sftp() as sftp_client:
				try: 
					sftp_client.stat(remotepath)
				except IOError: 
					print("[multibutler] File does not exist")
					return False
				else: 
					print("[multibutler] Downloading File...")
					sftp_client.get(remotepath, localpath)

					if os.path.isfile(localpath):
						return True
					else:
						return False
		except paramiko.SSHException:
			print("[multibutler] Connection Error")
			return False


	def download_file(self, localpath, tract, patch_s, filter, coadd='deepCoadd', filetype='calexp'): 

		# assemble path
		tail_path = super(self.__class__, self)._get_tail_path(tract=tract, patch_s=patch_s, filter=filter, coadd=coadd, filetype=filetype)
		remotepath = self.root_path + self.rerun_path + tail_path
		return self.download(remotepath, localpath)


	def download_psf(self, localpath, ra, dec, tract, patch_s, filter, coadd='deepCoadd', filetype='calexp'): 
		"""
		download psf 

		Params
		------
		localpath (string): path/file name to save to 
		ra (float)
		dec (float)
		tract (int): tract number
		patch_s (str): e.g., '7,3' patch_s number
		filter (str): e.g., 'HSC-G'
		coadd='deepCoadd'
		filetype='calexp'

		Return
		------
		status: True if successful, false if not.
		"""

		# assemble path
		tail_path = super(self.__class__, self)._get_tail_path(tract=tract, patch_s=patch_s, filter=filter, coadd=coadd, filetype=filetype)
		remotepath = self.root_path + self.rerun_path + tail_path

		dir_working = '/data/home/hscpipe/alsun/get_psf/'
		file_working = 'psf.fits'
		remotepsfpath = dir_working+file_working

		command1 = 'cd '+dir_working
		command2 = "ipython make_psf.py {remotepath} {file_working} {ra} {dec}".format(remotepath=remotepath, file_working=file_working, ra=str(ra), dec=str(dec))

		stdin, stdout, stderr = self.ssh_client.exec_command(command1+"; "+command2)

		try: 
			status_makepsf_insumire = (stdout.read().splitlines()[0] == 'True')
		except:
			status_makepsf_insumire = False

		if status_makepsf_insumire:
			status_download = self.download(remotepsfpath, localpath)
			stdin, stdout, stderr = self.ssh_client.exec_command('rm '+remotepsfpath)

			return status_download

		else:
			print("[multibutler] iaaButler.download_psf() make psf in sumire failed")
			print("[multibutler] please check if IAA_SUMIRE_CLUSTER_USERNAME is set to hscpipe")
			return False


	

	def _get_root_path(self):

		return '/array2/SSP/'



class onlineButler(parentButler):

	"""	

	Instruction
	-----------
	One needs to set sumire username and password as environment parameters

	$ export HSC_SSP_CAS_USERNAME
	$ read -s HSC_SSP_CAS_USERNAME
	$ export HSC_SSP_CAS_PASSWORD
	$ read -s HSC_SSP_CAS_PASSWORD

	"""
	def __init__(self, **kwargs):
		super(self.__class__, self).__init__(**kwargs)

		self.__username = getCrendential("HSC_SSP_CAS_USERNAME", cred_name = 'STARs username')
		self.__password = getCrendential("HSC_SSP_CAS_PASSWORD", cred_name = 'STARs password')

		self.root_path = self._get_root_path()


	def __enter__(self):
		""" nothing do be done"""
		pass


	def __exit__(self, exc_type, exc_value, traceback):
		""" nothing do be done"""
		pass


	def download_file(self, localpath, tract, patch_s, filter, coadd='deepCoadd', filetype='calexp'): 

		# assemble path
		tail_path = super(self.__class__, self)._get_tail_path(tract=tract, patch_s=patch_s, filter=filter, coadd=coadd, filetype=filetype)
		remotepath = self.root_path + self.rerun_path + tail_path

		# download
		rqst = requests.get(remotepath, auth=(self.__username, self.__password))
		if rqst.status_code == 200:

			with open(localpath, 'wb') as out:
				for bits in rqst.iter_content():
					out.write(bits)
			return True

		else:  
			print "image cannot be retrieved"
			return False


	def _get_root_path(self):
		return 'https://hscdata.mtk.nao.ac.jp/hsc_ssp/'



class sumireButler(parentButler):
	pass
