
class App
    constructor: () ->
        @data = new Data()
        #@pageNumber = new PageNumber(pages = 10, currentPage = 1)
        #@pageNumber2 = new PageNumber(pages = 20, currentPage = 5)
        @table = new Table(@data)
        #@filter = new Filter()
    view: ->
        [
            #try m("div", ["object1" ,@pageNumber.view()])
            #try m("div", ["object2" ,@pageNumber2.view()])
            #try m("div", ["object1" ,@pageNumber.view()])
            try header.view()
            try @filter.view()
            try @table.view()
            
        ]

class Data
    #model
    resp = {}
    meta=(response)->
        resp.headers = response.head.vars
        resp.genomes = response.results.bindings
        return
    constructor: () ->
        m.request(
            method: "POST",
            url: 'http://10.139.14.121:5000/mithril/meta',
            data: {}
            datatype: 'json'
            type: meta)
    #controller
    response: () ->
        return resp
    #view
    view: () ->
        return JSON.stringify(resp)

class Table
    ###
        Set $numbers to false to not display the left column in the table
        $pageY is window.pageYOffset (How far the user has scrolled down the page)
        $pageHeight is window.innerHeight (the area visable to the user.)
    ###
    numbers = true
    state = {pageY: 0, pageHeight: window.innerHeight}
    window.addEventListener("scroll", (e) ->
        state.pageY = Math.max(e.pageY || window.pageYOffset, 0)
        state.pageHeight = window.innerHeight
        m.redraw())
    constructor: (data) ->
        @data = data
    view: () ->
        data = @data.response()
        pageY = state.pageY
        begin = pageY / 46 | 0
        end = begin + (state.pageHeight /46 | 0 + 2)
        offset = pageY % 46
        m(".Occlusion", {style: {height: data.genomes.length * 46 + "px", position: "relative", top: -offset + "px"}}, [
            m("table", {style: {top: state.pageY + "px"}}, [
                m("tr", [
                    if numbers
                        m 'th' , "Redering: " + (begin * 1 + 1) + " to " + (end) + ". pageHeight: " + state.pageHeight#"#"
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



###
class Filter
    constructor: (data) ->
        @data = data
        @input_text= m.prop(2)
        @search_term= m.prop(5)
    view: -> [
        [
            m('input',{oninput: m.withAttr("value", @input_text), "value": @input_text()},[1]) 
            m('button',{onclick: m.withAttr("value", @search_term), "value": @input_text()},["Search"])
            
        ]
        m 'br'
        m('text',[@input_text()])
        m 'br'
        m('text',[@search_term()])
    ]
###

###
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
###