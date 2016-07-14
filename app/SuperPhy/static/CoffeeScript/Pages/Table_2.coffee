###
Table_2.coffee
    A trial table to test a different pagination style.
    URL: localhost/superphy/?/test

Args passed in:
    data: Data to populate the table
###
# coffeelint: disable=max_line_length

# Pages class: model for the table pages
class Pages

    constructor: (offset) ->
        this.offset = m.prop(offset) # This is the number of row displayed at a time
        this.currentpage = m.prop(1) # To display next to next and previous buttons
        this.numRows = m.prop(this.offset()-1) # A Class because ofm.prop. Sets the upper limit for number of rows to display at the time
        this.start = m.prop(0) # Row to begin display, adjusted in next()
    next: (max) ->
        if (max == false)
            return # The last row has been printed, do nothing
        else
            this.start(this.start() + this.offset()) # Started at 0, increments by 50
            this.currentpage(this.currentpage() + 1)
            this.numRows(this.numRows() + this.offset())
    prev: () ->    
        if (this.currentpage() == 1)
           return
        else
            this.numRows(this.numRows() - this.offset()) # Started at 50, decrements by 50 to change the offset in data stream
            this.start(this.start() - this.offset()) # Started at 0, decrements by 50 
            this.currentpage(this.currentpage() - 1)
    setpage: (num) ->
        this.currentpage = num



# Table class
class Table_2
    @controller: (args) ->
        @ViewPage = new Pages(50) # creating new page for the first page   
        
        # Sort function from Table.coffee which take two parameters: list and attribute
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
        this.rows = args.data().rows or m.prop([]) # The m.prop() here is an error handler, in case no rows are passed in 
        headers = args.data().headers
        #this.rows = debug_rows
        #headers = debug_headers
        
        len = @rows.length # this is the number of records in the data
        console.log("Number of records = #{len}, start = #{ctrl.ViewPage.start()}")
        lastRow = true
        pageNumber = Math.floor(len/ctrl.ViewPage.offset())+1 # This divides the number of records by the number of rows dispayed, and gets the number of pages
        console.log("#{pageNumber}")

        # This is to check if the maximum row has been printed, false if max reached
        if(ctrl.ViewPage.currentpage() >= pageNumber)
            lastRow = false
        # Overflow below was changed to auto from scroll.
        m('.', [   
            m("table", {"border": "1px solid black", "table-layout": "fixed", "overflow": "auto"}
                m("tr", [
                    for header in headers
                        m('th[data-sort-by=' + header + ']', ctrl.sort_table(list = this.rows) ,[header])
                        #m('th', {width: "150px", "text-align": "right"}, header) --this is code without sorting by headers
                ])
                
                for row, x in this.rows[ctrl.ViewPage.start() .. ctrl.ViewPage.numRows()] #when JSON.stringify(row).search(/Unknown/i) > -1
                    m('tr' 
                        for column in headers
                            m('td', {width: "110px"}, [row[column]])
                    )
            )
            if(ctrl.ViewPage.currentpage() != 1)
                m("button.btn btn-primary[type=button], left: 198px"
                    { onclick: -> ctrl.ViewPage.prev() } #Goes to the previous page method in Pages class
                        "Previous Page"
                )
            else
                m("button.btn"
                        onclick: -> ctrl.ViewPage.prev()
                        "Previous Page"
                )
            if(lastRow)
                m("button.btn btn-primary[type=button], right: 198px"
                        onclick: -> ctrl.ViewPage.next(lastRow)
                        "Next Page"
                )
            else
                m("button.btn"
                        onclick: -> ctrl.ViewPage.next(lastRow)
                        "Next Page"
                )
            m('p', ["Viewing page: #{ctrl.ViewPage.currentpage()} of #{pageNumber}"])
        ])

###
        table, th, td {
            "border": 1px solid black;
            "border-collapse": collapse;
        }
        th, td {
            padding: 15px;
        }
###