# GeneList component for gene search

class GeneList
    model: () =>
        @num_selected = 0

    constructor: () ->

        ## Things for the view (?)
        @state = {pageY: 0, pageHeight: window.innerHeight}
        window.addEventListener("scroll", (e) =>
            @state.pageY = Math.max(e.pageY || window.pageYOffset, 0)
            @state.pageHeight = window.innerHeight
            m.redraw()
        )


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
                        for header in headers.slice(0,1)
                            if row.visible() then \
                                m('td', {class: 'gene_table_item'}, [
                                    m('div', class: {'checkbox'}, [
                                        m('label', [
                                            m('input[type=checkbox]', {class: 'checkbox gene-table-checkbox gene-search-select', \
                                                                       checked: row.selected(), \
                                                                       onclick: m.withAttr("checked", row.selected)})
                                            [row[header]]
                                        ])
                                    ])
                                ])
                        for header in headers.slice(1)
                            if row.visible() then \
                                m('td', {class: 'gene_table_item'}, [
                                    m('label', [
                                        [row[header]]
                                    ])
                                ])
                    ])
            ])
        ])

