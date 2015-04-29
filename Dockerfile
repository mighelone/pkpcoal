FROM greole/python2.7

EXPOSE 9999

ADD . /root/pkp
WORKDIR /root/notebooks/pkp
RUN python setup.py install

CMD ip="$(ifconfig | grep -A 1 'eth0' | tail -1 | cut -d ':' -f 2 | cut -d ' ' -f 1)" && ipython notebook --ip $ip --port 9999
