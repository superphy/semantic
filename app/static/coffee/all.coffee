class App
    model: () =>
        @data ?= new Data()
        @table ?= new Table(@data)
        return
    controller: (options) =>
        @model()
        return
    view: () =>
        [
            #try m("div", ["object1" ,@pageNumber.view()])
            #try m("div", ["object2" ,@pageNumber2.view()])
            #try m("div", ["object1" ,@pageNumber.view()])
            header.view()
            @table.view()
        ]

class Data
    model: () =>
        meta=(response)=>
            @headers = response.head.vars
            @genomes = response.results.bindings
            return
        m.request(
            method: "POST",
            url: 'http://10.139.14.121:5000/mithril/meta',
            data: {}
            datatype: 'json'
            type: meta)
    #controller
    constructor: () ->
        @model()
    #view
    response: () =>
        return {
            headers: @headers
            genomes: @genomes
        }

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

sort_table= (list) ->
    { onclick: (e) ->
        item = e.target.getAttribute('data-sort-by')
        if item
            first = list[0]
            list.sort (a, b) ->
                if a[item]["value"] > b[item]["value"] then 1 else if a[item]["value"] < b[item]["value"] then -1 else 0
            if first == list[0]
                list.reverse()
        return
    }

class Home
    controller: ->

    #By default Coffeescript only returns the last item.
    #Wrap as an array to return all elements
    view: ->
        [
            header.view()
            m("div", {class:'container', id:'home-beta'}, [
                m("div", {class:'row'}, [
                    m("div", {class:'well center-block'}, [

                        m("p", {class:'text-center'}, [
                            m("span", {class:'text-info beta-release'}, "Beta Release")
                            ," Some features are still under development and may not be fully functional."
                        ])
                     ])
                ])
                m("div", {class:'row superphy-image'}, [
                    m("img", {src:'images/superphy_logo_with_title.png'}),
                    m("p", {class:'superphy-image'}, 'NEXT-LEVEL PHYLOGENETIC AND EPIDEMIOLOGICAL ANALYSIS OF PATHOGENS')
                ])
                m("div", {class:'row well'}, [
                    m("p", 'A user-friendly, integrated platform for the predictive genomic analyses of ',[
                        m("em", 'Escherichia coli')
                    ]),
                    m("p", 'Provides near real-time analyses of thousands of genomes using novel computational approaches.'),
                    m("p", 'Generates results that are understandable and useful to a wide community of users.')
                ])
                m("div", {class:'row'}, [
                    m("button", {class:'btn btn-danger btn-lg center-block'}, "INTRODUCTION")
                ])
                m("div", {class:'row text-center'}, [
                    m("p", 'For more information or to discuss possible collaboration, please contact:', [
                        m("p", [
                            m("ul", [
                                m("li", 'Dr. Vic Gannon: vic.gannon@phac-aspc.gc.ca')
                            ])
                        ])
                    ])
                ])
            ])
        ]

class Header

    #a list to be passed to the link generator function
    links: [
        {title: "Gnome", url: "/home"}
        {title: "Group Browse", url: "/gbrowse"}
        {title: "Group Analyses", url: "/groups"}
        {title: "VF and AMR", url: "/factors"}
        {title: "Meta", url: "/meta"}
    ]

    controller: ->

    view: ->
        m("div", {class:'container-fluid'}, [
            m("nav", {class:'navbar navbar-default navbar-fixed-top', role:'navigation'}, [
                m("div", {class:'navbar-header'}, [
                    m("a", {class:"navbar-brand", href:"/home", config:m.route}, "SuperPhy")
                ])
                m("ul", {class:'nav navbar-nav'}, [
                    m("li", m("a", {href: link.url, config:m.route}, link.title)) for link in @links
                    m("li", {class:"dropdown"}, [
                        m("a", {href:"", role:"button", class:"dropdown-toggle", 'data-toggle':"dropdown"}, "My Data", [
                            m("b", {class:"caret"})
                        ])
                    ])
                ])
            ])
        ])

header = new Header()
meta = new App

class Routes
    home = new Home
    m.route(document.body, "/", {
        "/": home
        "/home": home
        "/meta": meta
        "/gbrowse": meta
        "/groups": meta
    })