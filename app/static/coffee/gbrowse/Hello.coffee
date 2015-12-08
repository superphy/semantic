class Data
    resp = {}
    meta=(response)->
        resp.headers = response.head.vars
        resp.genomes = response.results.bindings
        return
    #controller
    constructor: () ->
        m.request(
            method: "POST",
            url: 'http://10.139.14.121:5000/mithril/meta',
            data: {}
            datatype: 'json'
            type: meta)
    data: () ->
        return resp

class PageNumber
    constructor: (pages = 1, currentPage = 1) ->
        @pages = pages
        @currentPage= m.prop(currentPage)
    next: ->
        @currentPage() * 1 + 1
    view: -> [
        [
            m('button',{onclick: m.withAttr("value", @currentPage), "value": x},[x]) 
            #if x % 10 == 0
            #    m("br")
        ] for x in [1 .. @pages]
        m('text',[@currentPage()])
    ]

class Table
    numbers = true
    constructor: (data) ->
        @data = data
    view: () ->
        data = @data.data()
        [
            m("table",[
                m("tr", [
                    if numbers
                        m 'th' , "#"
                    for header in data.headers
                        m('th[data-sort-by=' + header + ']',sort_table(data.genomes, header) ,[header]) 
                ])
                for binding, x in data.genomes
                    m("tr",[ 
                        if numbers
                            m 'td', x
                        for item in data.headers
                            try
                                m("td", binding[item]["value"])
                            catch
                                m("td", "")
                    ])
            ])
        ]

class App
    constructor: () ->
        @data = new Data()
        #@pageNumber = new PageNumber(pages = 10, currentPage = 1)
        #@pageNumber2 = new PageNumber(pages = 20, currentPage = 5)
        #@search = new Search()
        @table = new Table(@data)
    view: ->
        [
            try m("div", ["object1" ,@pageNumber.view()])
            try m("div", ["object2" ,@pageNumber2.view()])
            try m("div", ["object1" ,@pageNumber.view()])
            try header.view()
            try @search.view()
            try @table.view()
            
        ]

class Search
    constructor: () ->
        @currentPage= m.prop(2)
        @text= m.prop(5)
    next: ->
        @currentPage() * 1 + 1
    view: -> [
        [
            m('input',{oninput: m.withAttr("value", @currentPage), "value": @currentPage()},[1]) 
            m('button',{onclick: m.withAttr("value", @text), "value": @currentPage()},["Search"])
            
        ]
        m 'br'
        m('text',[@currentPage()])
        m 'br'
        m('text',[@text()])
    ]