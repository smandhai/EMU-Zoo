import numpy as np 
import argparse
import os

"HOW TO USE"

"""
Specify the txt file containing the DESTile names and the directory to store them in. See example below:

python DES_file_download.py /mnt/shared/home/user/DES_Tiles_EMU_ZOO.txt "/mnt/shared/des/DES_Tiles_new/"
"""
#Settings 
num =1 #2 # Can be values [1,2,3,4,5,6]
dr=1#Data Release

def main(args):
	#global DEStiles
	# '/Volumes/TARDIS/Work/EMUzoo/DES_tiles.txt'
	DESfile=str(args.infile[0])
	dataloc=str(args.dataloc[0])
	downloader=str(args.dataloc[0])

	rooturl='https://desdr-server.ncsa.illinois.edu/despublic/dr{}_tiles/'.format(dr)
	bands=['g',"i","r"] #choose subset of ['Y','g','i','r','z']
	#DEStiles=np.loadtxt(DESfile,dtype='str',delimiter='\n',skiprows=1)
	DEStiles=np.genfromtxt(DESfile,dtype='str',delimiter=' ',skip_header=1)
	print(DEStiles)
	for i in range(0,len(DEStiles)):
		tile=DEStiles[i]
		#tile = tile[0]
		print(f"Source: {tile}")
		for b in bands:
			url=rooturl+'{}/{}_r4575p01_{}.fits.fz'.format(tile,tile,b)
			#imagename='{}_r4575p01_{}.fits.fz'.format(tile,b)
			#curlstr="curl --header 'Host: desdr-server.ncsa.illinois.edu' --header 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/101.0.4951.64 Safari/537.36' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9' --header 'Accept-Language: en-US,en;q=0.9,en-GB;q=0.8' --header 'Referer: https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/{}' 'https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/{}/{}_r4575p01_{}.fits.fz' -L -o '{}{}_r4575p01_{}.fits.fz' ".format(tile,tile,tile,b,downloader,tile,b)
			if dr==1: #Use different naming ID if DR1 is used
				curlstr="curl -C - --fail --header 'Host: desdr-server.ncsa.illinois.edu' --header 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/101.0.4951.64 Safari/537.36' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9' --header 'Accept-Language: en-US,en;q=0.9,en-GB;q=0.8' --header 'Referer: https://desdr-server.ncsa.illinois.edu/despublic/dr{}_tiles/{}' 'https://desdr-server.ncsa.illinois.edu/despublic/dr{}_tiles/{}/{}_r2597p0{}_{}.fits.fz' -L -o '{}{}_r2590p0{}_{}.fits.fz' ".format(dr,tile,dr,tile,tile,num,b,downloader,tile,num,b)
			else:   #Use different naming format
				curlstr="curl -C - --fail --header 'Host: desdr-server.ncsa.illinois.edu' --header 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/101.0.4951.64 Safari/537.36' --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9' --header 'Accept-Language: en-US,en;q=0.9,en-GB;q=0.8' --header 'Referer: https://desdr-server.ncsa.illinois.edu/despublic/dr{}_tiles/{}' 'https://desdr-server.ncsa.illinois.edu/despublic/dr{}_tiles/{}/{}_r4575p0{}_{}.fits.fz' -L -o '{}{}_r4575p0{}_{}.fits.fz' ".format(dr,tile,dr,tile,tile,num,b,downloader,tile,num,b)
			os.system(curlstr)


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument('infile', type=str, nargs='+')
	ap.add_argument('dataloc', type=str, nargs='+')
	args = ap.parse_args()
	main(args)

