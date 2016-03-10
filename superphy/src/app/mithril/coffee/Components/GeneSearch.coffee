class GeneSearch
    constructor: () ->
        initGeneList = (gList, geneType, categories, tableElem, selElem, \countElem, catElem, autocElem, multi_select=true) =>
            for k,o of gList
                o.visible = true
                o.selected = false

            dObj = {
                type: geneType
                genes: gList
                categories: categories
                num_selected: 0
                sortField: 'name'
                sortAsc: true
                element: {
                    table: tableElem
                    select: selElem
                    count: countElem
                    category: catElem
                    autocomplete: autocElem
                }
                multi_select: multi_select
            }

        # FUNC appendGeneTable
        # Appends genes to form. Only attaches
        # genes where object.visible = true
        #
        # USAGE appendGeneList data_object
        #
        # RETURNS
        # boolean
        #  
        appendGeneTable = (d) =>
            table = d.element.table
            name = "#{d.type}-gene"
            tableElem = jQuery("<table />").appendTo(table)

            tableHtml = ''
            tableHtml += appendHeader(d)
            tableHtml += '<tbody>'
            tableHtml += appendGenes(d, sort(d.genes, d.sortField, d.sortAsc), d.genes, d.type, 'select')
            tableHtml += '</tbody>'

            tableElem.append(tableHtml)

            cboxes = table.find("input[name='#{name}']")
            cboxes.change( ->
                obj = $(@)
                geneId = obj.val()
                checked = obj.prop('checked')
                selectGene([geneId], checked, d)
            )
  
            updateCount(d)
  
            true

        @state = {pageY: 0, pageHeight: window.innerHeight}
        window.addEventListener("scroll", (e) =>
            @state.pageY = Math.max(e.pageY || window.pageYOffset, 0)
            @state.pageHeight = window.innerHeight
            m.redraw()
        )
        window.addEventListener("scroll", (e) =>
            @state.pageY = Math.max(e.pageY || window.pageYOffset, 0)
            @state.pageHeight = window.innerHeight
            m.redraw()
        )

    controller: () =>


    view: (data) =>
        rows = data.rows
        headers = data.headers
        pageY = @state.pageY
        begin = pageY / 46 | 0
        end = begin + (@state.pageHeight /46 | 0 + 10)
        offset = pageY % 46
        #console.log(pageY,begin,end,displaying)
        #height is the height of the entire table we want to display
        m(".Occlusion", {style: {height: rows.length * 46 + "px", position: "relative", top: -offset + "px"}}, [
            m("table", {style: {top: @state.pageY + "px"}}, [
                m("tr", [
                    for header in headers
                        m('th[data-sort-by=' + header + ']',events.sort_table(list = rows) ,[header])
                ])
                for row, x in rows[begin .. end] #when JSON.stringify(row).search(/Unknown/i) > -1
                    m('tr', [
                        for header in headers
                            m('td', [row[header]])
                    ])
            ])
        ])