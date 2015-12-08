
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
    state = {pageY: 0, pageHeight: window.innerHeight}
    window.addEventListener("scroll", (e) ->
        state.pageY = Math.max(e.pageY || window.pageYOffset, 0)
        state.pageHeight = window.innerHeight
        m.redraw())
    constructor: (data) ->
        @data = data
    view: () ->
        data = @data.data()
        pageY = state.pageY
        begin = pageY / 31 | 0
        end = begin + (state.pageHeight /31 | 0 + 2)
        offset = pageY % 31
        m(".Occlusion", {style: {height: data.genomes.length * 31 + "px", position: "relative", top: -offset + "px"}}, [
            m("table", {style: {top: state.pageY + "px"}}, [
                m("tr", [
                    if numbers
                        m 'th' , "Redering: " + (begin * 1 + 1) + " to " + (end * 1 + 1)#"#"
                    for header in data.headers
                        m('th[data-sort-by=' + header + ']',sort_table(data.genomes, header) ,[header]) 
                ])
                for binding, x in data.genomes[begin ... end]
                    m("tr",[ 
                        if numbers
                            m 'td', x * 1 + 1 + begin
                        for item in data.headers
                            try
                                m("td", binding[item]["value"])
                            catch
                                m("td", "")
                    ])
            ])
        ])


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