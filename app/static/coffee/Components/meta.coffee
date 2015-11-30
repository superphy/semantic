#Move headers and genomes to meta class
#Follow tutorial to get asych requests 

class Singleton
    @_instance: null
    @get:(arg)->

        @_instance or= new @()

class Meta_Data extends Singleton
    data = {}
    meta=(response)->
        data.headers = response.head.vars
        data.genomes = response.results.bindings
        return
    #controller
    constructor: ()->
        json = {data: {Genome_Uri:"F"}}
        m.request(
            method: "POST",
            url: 'http://10.139.14.121:5000/mithril/meta',
            data: json
            datatype: 'json'
            type: meta)
    views = {}
    views.table= () ->
        [
            m("div", {class:'container', id:'meta'},[
                m("table",[
                    m("tr", [
                        m("td", [item]) for item in data.headers
                    ])
                    for binding in data.genomes
                        m("tr",[ 
                            for item in data.headers
                                try
                                    m("td", binding[item]["value"])
                                catch
                                    m("td", "")
                        ])
                ])
            ])
        ]
    view: () ->
        [
            header.view()
            views.table()
        ]
