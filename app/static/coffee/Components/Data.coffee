class Data
    model: (json = {}) =>
        meta=(response)=>
            @headers = response.head.vars
            @rows = []
            for binding, x in response.results.bindings
                @rows[x] = {}
                for item in response.head.vars
                    try
                        @rows[x][item] = binding[item]['value']
                    catch
                        @rows[x][item] = ''
        m.request(
            method: "POST",
            url: 'http://' + location.hostname + ':5000/mithril/meta',
            data: json
            datatype: 'json'
            type: meta)
    #controller
    request: (json = {}) =>
        @model(json)
    constructor: (json = {}) ->
        @model(json)
    #view
    response: () =>
        return {
            rows: @rows
            headers: @headers
        }