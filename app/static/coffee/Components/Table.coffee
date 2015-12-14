class Table
    model: (data) =>
        ###
            Set $numbers to false to not display the left column in the table
            $pageY is window.pageYOffset (How far the user has scrolled down the page)
            $pageHeight is window.innerHeight (the area visable to the user.)
        ###
        @numbers = true
        @state = {pageY: 0, pageHeight: window.innerHeight}
        @data = data

    constructor: (data) ->
        @model(data)
        window.addEventListener("scroll", (e) =>
            @state.pageY = Math.max(e.pageY || window.pageYOffset, 0)
            @state.pageHeight = window.innerHeight
            m.redraw())
        window.addEventListener("resize", (e) =>
            @state.pageY = Math.max(e.pageY || window.pageYOffset, 0)
            @state.pageHeight = window.innerHeight
            m.redraw())
        
    view: () ->
        data = @data.response()
        pageY = @state.pageY
        begin = pageY / 46 | 0
        end = begin + (@state.pageHeight /46 | 0 + 2)
        offset = pageY % 46
        m(".Occlusion", {style: {height: data.genomes.length * 46 + "px", position: "relative", top: -offset + "px"}}, [
            m("table", {style: {top: @state.pageY + "px"}}, [
                m("tr", [
                    if @numbers
                        m 'th' , "Redering: " + (begin * 1 + 1) + " to " + (end) + ". pageHeight: " + @state.pageHeight#"#"
                    for header in data.headers
                        m('th[data-sort-by=' + header + ']',events.sort_table(list = data.genomes) ,[header]) 
                ])
                for binding, x in data.genomes[begin ... end]
                    m("tr",[ 
                        if @numbers
                            m 'td', x * 1 + 1 + begin
                        for item in data.headers
                            try
                                m("td", binding[item]["value"])
                            catch
                                m("td", "")
                    ])
            ])
        ])