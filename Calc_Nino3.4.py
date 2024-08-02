def StandardGrid(data,lon_name='longitude',lat_name='latitude'):
	"""
		Make sure longitude is in [0,360] and latitude sorted
		 from lowest to highest.
		 Assumes an xr.DataArray or xr.Dataset.

		INPUTS:
			data:	   xarray.DataArray or xarray.Dataset to regrid
			lon_name:  name of longitude dimension. Set to None if nothing should be done.
			lat_name:  name of latitude dimension. Set to None if nothing should be done.
		OUTPUTS:
			data:	  xarray.DataArray with latitude from lowest to highest and
					   longitude between 0 and 360 degrees.
	"""
	if lat_name is not None and lat_name in data.coords:
		if data[lat_name][0] > data[lat_name][-1]:
			data = data.sortby(lat_name)
	if lon_name is not None and lon_name in data.coords and data[lon_name].min() < 0:
		data = data.assign_coords({lon_name : (data[lon_name]+360)%360})
		return data.sortby(lon_name)
	else:
		return data
    

def Nino(sst, lon='lon', lat='lat', time='time', avg=5, nino='3.4'):
	"""
		Produce ENSO index timeseries from SST according to Technical Notes
		 guidance from UCAR: https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni

		INPUTS:
			sst:  xarray.DataArray which will be averaged over Nino domains
			lon:  name of longitude dimension. Has to be in [0,360].
			lat:  name of latitude dimension. Has to be increasing.
			time: name of time dimension.
			avg:  size of rolling window for rolling time average.
			nino: which Nino index to compute. Choices are
					'1+2','3','4','3.4','oni','tni'

		OUTPUTS:
			sst: spatially averaged over respective Nino index domain
				  note that no running means are performed.
	"""
	ninos = {
		'1+2' : {lon:slice(270,280),lat:slice(-10,0)},
		'3'   : {lon:slice(210,270),lat:slice(-5,5)},
		'4'   : {lon:slice(160,210),lat:slice(-5,5)},
		'3.4' : {lon:slice(190,240),lat:slice(-5,5)},
		'oni' : {lon:slice(190,240),lat:slice(-5,5)},
	}
	possible_ninos = list(ninos.keys())+['tni']
	if nino not in possible_ninos:
		raise ValueError('Nino type {0} not recognised. Possible choices are {1}'.format(nino,', '.join(possible_ninos)))
	lon_name = None
	lat_name = None
	if sst[lon].min() < 0 or sst[lon].max() <= 180:
		lon_name = lon
	if sst[lat][0] > sst[lat][-1]:
		lat_name = lat
	if lon_name is not None or lat_name is not None:
		print('WARNING: re-arranging SST to be in domain [0,360] x [-90,90]')
		sst = StandardGrid(sst,lon_name,lat_name)

	def NinoAvg(sst,nino,time,avg):
		ssta = sst.sel(ninos[nino]).mean(dim=[lon,lat])
		sstc = ssta.groupby('.'.join([time,'month'])).mean(dim=time)
		ssta = ssta.groupby('.'.join([time,'month'])) - sstc
		if avg is not None:
			ssta = ssta.rolling({time:avg}).mean()
		return ssta/ssta.std(dim=time)

	if nino == 'tni':
		n12 = NinoAvg(sst,'1+2',time,None)
		n4  = NinoAvg(sst,'4',time,None)
		tni = (n12-n4).rolling({time:avg}).mean()
		return tni/tni.std(dim=time)
	else:
		return NinoAvg(sst,nino,time,avg)
    