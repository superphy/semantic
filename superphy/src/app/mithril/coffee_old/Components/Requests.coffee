class Request
    #ex: model: @data ?= new Request(url="data/meta", type='sparql')
    #ex: view : @data.view()
    model: () =>
        m.request(
            method: @method,
            url: "http://#{location.hostname}:5000/#{@url}",
            data: @data
            datatype: @datatype
            type: (response)=>
                switch @type
                    when 'array'
                        #expect data as {data: [ {headers : values}]}
                        #expect missing keys to be present with '' values
                        @rows = response.data
                        @headers = []
                        for key of @rows[0]
                            @headers.push(key)
                    when 'sparql'
                        #expect data as {head:{vars:[], results:{bindings:[]}}
                        @headers = response.head.vars
                        @rows = []
                        for binding, x in response.results.bindings
                            @rows[x] = {}
                            for item in response.head.vars
                                try
                                    @rows[x][item] = binding[item]['value']
                                catch
                                    @rows[x][item] = ''
        )
    constructor: (url, type='array', method='GET', data={}, datatype='json') ->
        if url.substring(0, 1) == '/'
            url = url.substring(1)
        @url = url
        @type ?= type.toLowerCase()
        @method ?= method.toUpperCase()
        @data ?= data
        @datatype ?= datatype
        @model()

    table_view: () =>
        m("table", [
            m("tr", [
                for header in @headers
                    m('th', [header])
            ])
            for row, x in @rows[0 .. 10]
                m('tr', [
                    for header in @headers
                        m('td', [row[header]])
                ])
        ])
        
    view: () =>
        @table_view()

class RequestMeta extends Request
    #ex: model: @data ?= new RequestMeta()
    #ex: view: @data.view()
    constructor: () ->
        @url ?= "data/meta"
        @type ?= "sparql"
        super(@url)