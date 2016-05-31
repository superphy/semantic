# coffeelint: disable=max_line_length

MetaData = (args) ->
    args = args || {}
    #at this point the data is static. At a later date we can add in
    #expire times, or a query to the server to ask if it should be updated.
    try
        @data ?= m.prop(JSON.parse(localStorage.getItem("MetaData")))
    catch
        @data = m.prop(null)
    if @data() is null or args.reset is true
        @data = m.request(
            method: "POST",
            url: "http://#{location.hostname}:5000/data/meta"
            data: {}
            datatype: 'json',
            type: (response) ->
                data = {}
                data.headers = response.head.vars
                data.rows = []
                for binding, x in response.results.bindings
                    data.rows[x] = {}
                    for item in response.head.vars
                        try
                            data.rows[x][item] = binding[item]['value']
                        catch
                            data.rows[x][item] = ''
                data.date = new Date()
                return data
        )
    window.onbeforeunload = ->
        localStorage.setItem("MetaData", JSON.stringify(@data()))
    return @data