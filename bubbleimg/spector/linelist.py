# linelist.py
# ALS 2017/06/14

""" include line list for spectral line related operations """


# ======= for continuum subtraction ===========
narrowline = ['OII', 'OIII', 'OI', 'NII', 'NeV', 'SII', 'Ha', 'Hb', 'Hg', 'Hd', 'ArIII', 'NeIII', 'Dip', 'Mg', 'Na', 'He', 'HeI', 'HeII', 'FeVII']

broadline = ['Ha', 'Hb', 'Hg', 'Hd']

#======== for inferring line flux =============
# strongline as defined by lines that are brighter than 1/50 of OIII 5008
strongline = ['MgII2800', 'OII3730', 'NeIII3870', 'NeIII3969', 'Hg', 'Hb', 'OIII4960', 'OIII5008', 'OI6302', 'NII6550', 'Ha', 'NII6585', 'SII6718', 'SII6733', ]
# 'OII3726', excluded as it is very close to OII3730
# 'OIII4364' marginal. almost fainter than 2% of [OIII]
# Caveat: NeV 3427 is not a very faint line but is not measured by SDSS


# mapping from my line tag to sdss LINENAME, not all the sdss lines are included. 
# H_epsilon of SDSS is not included because I suspect sdss uses a wrong wavelength of 3890 while the true value is 3970. 
sdssLINENAME = {
  'MgII2800' : 'Mg_II 2799', 
  'OII3726' : '[O_II] 3725', 
  'OII3730' : '[O_II] 3727', 
  'NeIII3870' : '[Ne_III] 3868', 
  'NeIII3969' : '[Ne_III] 3970', 
  'Hg' : 'H_gamma', 
  'OIII4364' : '[O_III] 4363', 
  'Hb' : 'H_beta', 
  'OIII4960' : '[O_III] 4959', 
  'OIII5008' : '[O_III] 5007', 
  'OI6302' : '[O_I] 6300', 
  'OI6366' : '[O_I] 6363', 
  'NII6550' : '[N_II] 6548', 
  'Ha' : 'H_alpha', 
  'NII6585' : '[N_II] 6583', 
  'SII6718' : '[S_II] 6716', 
  'SII6733' : '[S_II] 6730', 
  'ArIII7138' : '[Ar_III] 7135', 
  }