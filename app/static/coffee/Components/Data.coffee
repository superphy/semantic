class Data
    model: (json = {}) =>
        meta=(response)=>
            @headers = response.head.vars
            @genomes = response.results.bindings
            return
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
            headers: @headers
            genomes: @genomes
        }
