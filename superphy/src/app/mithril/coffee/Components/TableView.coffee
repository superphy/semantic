###
CLASS TableView

Component class for table views

Args passed in:
    data: Data to populate the table
    checkbox: boolean for checkbox option (for VF & AMR)

To fix: make the VF & AMR checkbox portion look nicer?
###

class TableView
    state = {pageY: 0, pageHeight: window.innerHeight}
    window.addEventListener("scroll", (e) ->
        state.pageY = Math.max(e.pageY || window.pageYOffset, 0)
        state.pageHeight = window.innerHeight
        m.redraw()
    )
    @controller: (args) ->
        @sort_table = (list, attribute = 'data-sort-by') ->
            { onclick: (e) ->
                item = e.target.getAttribute(attribute)
                if item
                    first = list[0]
                    list.sort (a, b) ->
                        if isNaN(parseFloat(a[item] * 1))
                            if isNaN(parseFloat(b[item] * 1))
                                if a[item] > b[item] then 1
                                else if b[item] > a[item] then -1
                                else 0
                            else -1
                        else if isNaN(parseFloat(b[item] * 1)) then 1
                        else if a[item] * 1 < b[item] * 1 then 1
                        else if b[item] * 1 < a[item] * 1 then -1
                        else 0
                    if first == list[0]
                        list.reverse()
                return
            }
        return @
    
    @view: (ctrl, args) ->
        rows = args.data.rows
        headers = args.data.headers
        pageY = state.pageY
        begin = pageY / 60 | 0
        end = begin + (state.pageHeight /60 | 0 + 10)
        offset = pageY % 60
        m(".Occlusion", {style: {height: args.data.rows.length * 46 + "px", position: "relative", top: -offset + "px"}}, [
            m("table", {style: {top: state.pageY + "px"}}, [
                m("tr", [
                    for header in args.data.headers
                        m('th[data-sort-by=' + header + ']', ctrl.sort_table(list = args.data.rows) ,[header])
                ])
                for row, x in args.data.rows[begin .. end] #when JSON.stringify(row).search(/Unknown/i) > -1
                    m('tr', [
                        if args.checkbox and row.visible()
                            for header in headers
                                ## For VF and AMR
                                if header is "Gene"
                                    m('td', {class: 'gene_table_item'}, [
                                        m('.checkbox', [
                                            m('label', [
                                                m('input[type=checkbox]', 
                                                  {class: 'checkbox gene-table-checkbox gene-search-select', \
                                                   checked: row.selected(), \
                                                   onclick: m.withAttr("checked", row.selected)})
                                                row[header]
                                            ])
                                        ])
                                    ])
                                else
                                    m('td', {class: 'gene_table_item'}, [
                                        m('label', [
                                            row[header]
                                        ])
                                    ])
                        else if not row.visible()
                            null
                        else
                            for header in args.data.headers
                                m('td', [row[header]])
                    ])
            ])
        ])