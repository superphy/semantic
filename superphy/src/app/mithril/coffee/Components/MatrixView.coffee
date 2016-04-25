class MatrixView
    constructor: () ->
        return @

    init: (searchResults, parentElem, elID) =>
        ## search results is an object
        # console.log("object", searchResults)
        # console.log("object keys", Object.keys(searchResults))
        genomes = Object.keys(searchResults)
        genes = Object.keys(searchResults[genomes[0]])

        view = @

        (el, isInitialized, ctx) ->
            if not isInitialized
                view._create_matrix(genomes, genes, searchResults, parentElem, elID)

    _create_matrix: (genomes, genes, searchResults, parentElem, elID) =>
        console.log("Creating matrix...")
        self = @

        @cellWidth = 20
        @margin = {top: 150, right: 0, bottom: 0, left:250}
        @height = genomes.length * @cellWidth
        @width = genes.length * @cellWidth
        dim = {
            w: @width + @margin.right + @margin.left
            h: @height + @margin.top + @margin.bottom
        }

        ## parent elements
        @parentElem = parentElem
        @elID = elID


        # sets @genomeNodes, @geneNodes, @matrix, @lengenomes, @lengenes
        @_compute_matrix(genomes, genes, searchResults)


        # Precompute ordering for genes which are stable
        # genomes can be removed through filtering
        @geneOrders = {
            name: d3.range(@lengenes).sort (a, b) =>
                d3.ascending(@geneNodes[a].name, @geneNodes[b].name)

            count: d3.range(@lengenes).sort (a, b) =>
                @geneNodes[b].count - @geneNodes[a].count
        }
        @geneOrders['group'] = @geneOrders['count']
        @orderType = 'name' #default is alphabetical

        # Set color opacity based on num of alleles
        # Above 4 alleles full saturation is reached
        @z = d3.scale.linear().domain([0, 4]).clamp(true)
        @x = d3.scale.ordinal().rangeBands([0, @width])

        @cssClass = 'superphy-matrix'

        ## Create sort drop down
        dropdownDiv = jQuery('<div class="matrixSort"><span>Order:</span> </div>').appendTo("##{@elID}")
        dropdown = jQuery('<select name="matrix-sort">' +
                '<option value="name" selected="selected"> by Name</option>' +
                '<option value="count"> by Frequency</option>' +
                '<option value="group"> by Group</option>' +
                '</select>').appendTo(dropdownDiv)

        ## view number (from Matt's previous work)
        num = 0 # @elNum - 1
        dropdown.change( ->
            sortType = @.value
            console.log("sortType", sortType)
            self.viewAction(num, ['matrix-sort', sortType])
        )

        wrap = d3.select("##{@elID}")
            .append("div")
            .attr("class", "matrix-container")
            .append("svg")
            .attr("width", dim.w)
            .attr("height", dim.h)

        ## Parent element for attaching matrix elements
        @canvas = wrap.append("g")
            .attr("transform", "translate(" + @margin.left + "," + @margin.top + ")")

        ## Background
        @canvas.append("rect")
            .attr("class", "matrixBackground")
            .attr("width", @width)
            .attr("height", @height)
            .style("fill", "#f2f2f2")

        type: "matrix"
        elName: "genome_matrix"
        duration: 500

        @formatCount = d3.format(",.0f")

        ## initial load
        self.viewAction(num, ['matrix-sort', 'name'])

    _compute_matrix: (genomes, genes, results) =>
        @lengenomes = genomes.length
        @lengenes = genes.length

        # Genome nodes
        @matrix = []

        @genomeNodes = []
        i = 0
        for g in genomes
            gObj = {
                id: i
                name: g
                count: 0
            }
            @genomeNodes.push(gObj)
            @matrix[i] = d3.range(@lengenes).map (j) ->
                {x: j, y: i, z: 0, i:null}
            i++

        @geneNodes = []
        i = 0
        for g in genes
            gObj = {
                id: i
                name: g
                count: 0
            }
            @geneNodes.push(gObj)
            i++

        i = 0
        for g in @genomeNodes
            for n in @geneNodes
                numAlleles = results[g.name][n.name]

                g.count += numAlleles
                n.count += numAlleles

                @matrix[g.id][n.id].z = numAlleles
                @matrix[g.id][n.id].i = i

                i++

    viewAction: (genomes, argArray) ->
        evt = argArray.shift()

        if evt is 'matrix-sort'
            @orderType = argArray[0]
            throw new SuperphyError "Unrecognized order type: #{@orderType} in MatrixView viewAction method." unless @orderType in Object.keys(@geneOrders)
            @_update(genomes)
          
        else
          throw new SuperphyError "Unrecognized event type: #{evt} in MatrixView viewAction method."


    _update: (genomes) ->
        t1 = new Date()

        @_sync(genomes)

        @height = @cellWidth * @lengenomes

        @y = d3.scale.ordinal().rangeBands([0, @height])
        @y.domain(@genomeOrders[@orderType])
        @x.domain(@geneOrders[@orderType])

        @canvas.selectAll(".matrixBackground")
            .attr("height", @height)

        # Attach genome rows
        svgGenomes = @canvas.selectAll("g.matrixrow")
            .data(@genomeNodes, (d) -> d.id)

        svgGenomes
            .attr("class", (d) => @_classList(d))
            .select("text.matrixlabel")
            .text((d) -> d.name)

        svgGenomes
            .selectAll("g.matrixcell title")
            .text((d) -> d.name + "title")

        # Insert new rows at origin
        that = @
        newRows = svgGenomes.enter().append("g")
            .attr("class", (d) => @_classList(d))
            .attr("transform", (d, i) -> "translate(0,0)")
            .each((d) -> that._row(@, that.matrix[d.id], that.x, that.y, that.z))

        newRows.append("line")
            .attr("x2", @width)
            .style("fill", "white")

        newRows.append("text")
            .attr("class","matrixlabel")
            .attr("x", -6)
            .attr("y", @y.rangeBand() / 2)
            .attr("dy", ".32em")
            .attr("text-anchor", "end")
            .text((d) -> d.name)
            .on("click", (d) ->
                jQuery('#dialog-matrix-row-select')
                    .data('row-id', d.genome)
                    .dialog('open')
            )

        svgGenes = @canvas.selectAll("g.matrixcolumn")
            .data(@geneNodes, (d) -> d.id)

        svgGenes.selectAll("line")
            .attr("x1", -@height)
            .style("color", "white")

        # Insert new columns at origin
        newCols = svgGenes.enter().append("g")
            .attr("class", "matrixcolumn")
            .attr("transform", (d, i) -> "translate(" + 0 + ")rotate(-90)" )

        newCols.append("line")
            .attr("x1", -@height)
            .style("color", "white")

        newCols.append("text")
            .attr("class","matrixlabel")
            .attr("x", 6)
            .attr("y", @y.rangeBand() / 2)
            .attr("dy", ".32em")
            .attr("text-anchor", "start")
            .text((d) -> d.name)
            .on("click", (d) ->
                jQuery('#dialog-matrix-col-select')
                    .data('col-id', d.gene)
                    .dialog('open')
            )

        # Move new and existing rows to final position
        @_assumePositions()

        genomesExit = svgGenomes.exit().transition()
            .duration(@duration)
            .attr("transform", (d) -> "translate(0,0)")
            .remove()

    _sync: (genomes) ->

        @currNodes = []
        @currN = 0

        for n in @genomeNodes
            n.viewname = n.name
            n.assignedGroup = 0

            n.index = @currN
            @currNodes.push(n)
            @currN++

            i = 0
            for c in @matrix[n.id]
                c.title = "genome title"
                i++

        @genomeOrders = {
            name: d3.range(@currN).sort (a, b) =>
                d3.ascending(@currNodes[a].viewname, @currNodes[b].viewname)
                
            count: d3.range(@currN).sort (a, b) =>
                @currNodes[b].count - @currNodes[a].count
              
            group: d3.range(@currN).sort (a, b) =>
                gdiff = @currNodes[b].assignedGroup - @currNodes[a].assignedGroup
                if gdiff == 0
                    return @currNodes[b].count - @currNodes[a].count
                else
                    return gdiff
        }
        

    _row: (svgRow, rowData, x, y, z) ->
        #Grab row cells
        svgCells = d3.select(svgRow).selectAll("matrixcell")
            .data(rowData, (d) -> d.i)

        # Insert new cells
        num = 0 #@elNum-1

        newCells = svgCells.enter().append("g")
            .attr("class", "matrixcell")
            .attr("transform", (d, i) =>
                "translate(" + @x(d.x) + ",0)")
            .on("mouseover", (p) => @_mouseover(p))
            .on("mouseout", @_mouseout)

        newCells.append("rect")
            .attr("x", 0)
            .attr("width", x.rangeBand())
            .attr("height", y.rangeBand())
            .style("fill-opacity", (d) -> z(d.z))

        newCells.append("text")
            .attr("dy", ".32em")
            .attr("y", x.rangeBand() / 2)
            .attr("x", x.rangeBand() / 2)
            .attr("text-anchor", "middle")
            .text( (d) =>
                if d.z > 0
                    @formatCount(d.z)
                else
                    ''
            )
        true

    _assumePositions: ->
    
        that = @
        @duration = 250 ## default value
        transit = @canvas.transition().duration(@duration)
        #transit = @canvas.transition().
        
        transit.selectAll(".matrixrow")
            #.delay((d, i) -> that.y(d.index) * 4)
            .attr("transform", (d, i) ->
                "translate(0," + that.y(d.index) + ")")
            .selectAll(".matrixcell")
            #.delay((d) -> that.x(d.x) * 4)
            .attr("x", (d) -> that.x(d.x))
            .attr("transform", (d, i) ->
                "translate(" + that.x(d.x) + ",0)")

        transit.selectAll(".matrixcolumn")
            #.delay((d, i) -> that.x(i) * 4)
            .attr("transform", (d, i) ->
                "translate(" + that.x(i) + ")rotate(-90)")
            
        true

    _mouseover: (p) ->
        d3.selectAll(".matrixrow text")
            .classed("matrixActive", (d, i) -> d.index == p.y)
        d3.selectAll(".matrixcolumn text")
            .classed("matrixActive", (d, i) -> i == p.x)

    _mouseout: ->
        d3.selectAll("text").classed("matrixActive", false)


    _classList: (d) ->
        clsList = ['matrixrow']
        clsList.push("selectedRow") if d.selected
        clsList.push("groupedRow#{d.assignedGroup}") if d.assignedGroup?

        clsList.join(' ')