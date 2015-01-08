# -*- coding: utf-8 -*-
"""
"""
import os
import redis
import urlparse
from werkzeug.wrappers import Request, Response
from werkzeug.routing import Map, Rule
from werkzeug.exceptions import HTTPException, NotFound
from werkzeug.wsgi import SharedDataMiddleware
from werkzeug.utils import redirect

from jinja2 import Environment, FileSystemLoader

from PKP import PKP

def base36_encode(number):
    assert number >= 0, 'positive integer required'
    if number == 0:
        return '0'
    base36 = []
    while number != 0:
        number, i = divmod(number, 36)
        base36.append('0123456789abcdefghijklmnopqrstuvwxyz'[i])
    return ''.join(reversed(base36))

def get_hostname(url):
    return urlparse.urlparse(url).netloc

class Shortly(object):

    def __init__(self, config):
        self.redis = redis.Redis(config['redis_host'], config['redis_port'])
        template_path = os.path.join(os.path.dirname(__file__), 'templates')
        self.jinja_env = Environment(loader=FileSystemLoader(template_path),
                                     autoescape=True)
        self.jinja_env.filters['hostname'] = get_hostname

        self.url_map = Map([
            Rule('/', endpoint='inputs'),
            Rule('/preProcDict', endpoint='preProcDict'),
            Rule('/preProcList', endpoint='preProcList'),
            Rule('/preProcTSV', endpoint='preProcTSV'),
            Rule('/model', endpoint='short_link_details')
        ])

    def on_inputs(self, request):
        import os
        error = None
        url = ''
        if request.method == 'POST':
            url = request.form
            input_dict = self.transform_dict(dict(url))
            print input_dict
            generator = PKP.Generate(input_dict)
            self.res = generator.executeSolver()[0] #FIXME enable multiple runs
            res = str(self.res._tsv)
            target =  os.getcwd() + "/static/res.tsv"
            with open(target,'w') as f:
                l = res.replace('(ms)','').replace(' ','\t')
                f.write(l)
        return self.render_template('inputs.html', error=error, url=url)

    def transform_dict(self, raw_dict):
        return {'Coal': 
        {'Ultimate Analysis'  : {
            'Carbon'              :float(raw_dict['carbon'][0]),
            'Hydrogen'            :float(raw_dict['hydrogen'][0]),
            'Oxygen'              :float(raw_dict['oxygen'][0]),
            'Nitrogen'            :float(raw_dict['nitrogen'][0])},
         'Proximate Analysis' : {
            'Volatile Matter'     :float(raw_dict['volatile_matter'][0]),
            'Moisture'            :float(raw_dict['moisture'][0]),
            'Fixed Carbon'        :float(raw_dict['fixed_carbon'][0]),
            'Ash'                 :float(raw_dict['ash'][0])}},
        'CPD' : {
            'active': True,
            'deltaT':   1e-5,
        },
        'OperatingConditions':{
            'pressure' : 1.0,
            'runs'     : 1,
            'run1'     : [ [ 0, 1000], [ 0.5,1500], [ 1.0, 2000] ]
            }
        }

    def on_preProcDict(self, request):
        #res = PKP.generate(json_string=os.getcwd() + '/../../inputs/')[0].data
        res = str(self.res._dict).replace('\'','')
        return self.render_template('json.html', error=None, title= res)

    def on_preProcList(self, request):
        #res = PKP.generate(json_string=os.getcwd() + '/../../inputs/')[0].data
        res = str(self.res._list).replace('\'','')
        return self.render_template('json.html', error=None, title= res)

    def on_preProcTSV(self, request):
        #res = PKP.generate(json_string=os.getcwd() + '/../../inputs/')[0].data
        res = str(self.res._tsv)
        target =  os.getcwd() + "static/res.tsv"
        
        return self.render_template('json.html', error=None, title= res)

    def error_404(self):
        response = self.render_template('404.html')
        response.status_code = 404
        return response


    def render_template(self, template_name, **context):
        t = self.jinja_env.get_template(template_name)
        return Response(t.render(context), mimetype='text/html')

    def dispatch_request(self, request):
        adapter = self.url_map.bind_to_environ(request.environ)
        try:
            endpoint, values = adapter.match()
            return getattr(self, 'on_' + endpoint)(request, **values)
        except NotFound, e:
            return self.error_404()
        except HTTPException, e:
            return e

    def wsgi_app(self, environ, start_response):
        request = Request(environ)
        response = self.dispatch_request(request)
        return response(environ, start_response)

    def __call__(self, environ, start_response):
        return self.wsgi_app(environ, start_response)


def create_app(redis_host='localhost', redis_port=6379, with_static=True):
    app = Shortly({
        'redis_host':       redis_host,
        'redis_port':       redis_port
    })
    if with_static:
        app.wsgi_app = SharedDataMiddleware(app.wsgi_app, {
            '/static':  os.path.join(os.path.dirname(__file__), 'static')
        })
    return app


if __name__ == '__main__':
    from werkzeug.serving import run_simple
    app = create_app()
    run_simple('127.0.0.1', 5001, app, use_debugger=True, use_reloader=True)
