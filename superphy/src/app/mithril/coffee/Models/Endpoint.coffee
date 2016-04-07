# coffeelint: disable=max_line_length

getEndpoint = (url) ->
    try
        response = m.prop(JSON.parse(localStorage.getItem(url)))
        diff = (new Date(response().date).getTime() - new Date().getTime())
        #console.log(diff)
        if diff < 0
            throw diff
    catch
        localStorage.setItem(url, null)
        response = m.prop(null)
    if response() is null
        response = m.request(
            method: "POST",
            url: "http://#{location.hostname}:5000/#{url}"
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
                data.date = response.date
                localStorage.setItem(url, JSON.stringify(data))
                return data
        )
    return response

class RawEndpoint
    @controller: (args) ->
        url = args.url
        @data = getEndpoint(args.url)
        return @
    @view: (ctrl)->
        m('.', JSON.stringify(ctrl.data()))