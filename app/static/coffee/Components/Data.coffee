class Data
    model: (json = {}) =>
        noID=(response)=>
            @_headers = response.head.vars
            @_rows = []
            for binding, x in response.results.bindings
                @_rows[x] = {}
                for item in response.head.vars
                    try
                        @_rows[x][item] = binding[item]['value']
                    catch
                        @_rows[x][item] = ''
            @search('')
        ID=(response)=>
            @_headers = ["id"].concat(response.head.vars)
            @_rows = []
            for binding, x in response.results.bindings
                @_rows[x] = {}
                @_rows[x]["id"] = x
                for item in response.head.vars
                    try
                        @_rows[x][item] = binding[item]['value']
                    catch
                        @_rows[x][item] = ''
            @search('')
        m.request(
            method: "POST",
            url: 'http://' + location.hostname + ':5000/mithril/meta',
            data: json
            datatype: 'json'
            type: ID)
    #controller
    request: (json = {}) =>
        @model(json)
    constructor: (json = {}) ->
        @model(json)
    #view
    search: (searchterm) =>
        searchterm = searchterm.toLowerCase()
        @rows = []
        @rows.push(row) for row in @_rows when JSON.stringify(row).toLowerCase().search(searchterm) > -1
        console.log(searchterm)
        @headers = @_headers
    response: () =>
        return {
            rows: @rows
            headers: @headers
        }