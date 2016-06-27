
#This does not follow coding conventions. A refactor is in order.

class HistogramView
    constructor: () ->
        return @

    init: (searchResults, parentElem, elID) =>
        ## search results is an object
        genomes = Object.keys(searchResults)
        genes = Object.keys(searchResults[genomes[0]])

        view = @

        (el, isInitialized, ctx) ->
            #console.log(ctx)
            if not isInitialized
                view._create_histogram(searchResults, parentElem, elID)


    _create_histogram: (searchResults, parentElem, elID) =>
        ## Create empty histogram
        margin = {top: 40, right: 30, bottom: 40, left: 30}
        @width = 300 - margin.left - margin.right
        @height = 250 - margin.top - margin.bottom

        ## parent elements
        @parentElem = parentElem
        @elID = elID
        @cssClass = "matrix_histogram"

        bins = [
            {'val': 0, 'key': '0'}
            {'val': 1, 'key': '1'}
            {'val': 2, 'key': '2'}
            {'val': 3, 'key': '3'}
            {'val': 4, 'key': '4'}
            {'val': 5, 'key': '>=5'}
        ]

        @x = d3.scale.ordinal()
            .domain(bins.map((d) -> d.val))
            .rangeRoundBands([0, @width], .05)

        @x2 = d3.scale.ordinal()
            .domain(bins.map((d) -> d.key))
            .rangeRoundBands([0, @width], .05)

        @xAxis = d3.svg.axis()
            .scale(@x2)
            .orient("bottom")

        @histogram = d3.layout.histogram()
            .bins([0,1,2,3,4,5,6])

      
        @canvas = d3.select("##{@elID}").append("svg")
            .attr("width", @width + margin.left + margin.right)
            .attr("height", @height + margin.top + margin.bottom)
            .append("g")
            .attr("transform", "translate(#{margin.left},#{margin.top})")
      
        @formatCount = d3.format(",.0f")
    
        @canvas.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(0," + @height + ")")
            .call(@xAxis)
            .append("text")
            .attr("dy", ".75em")
            .attr("y", 23)
            .attr("x", @width / 2)
            .attr("text-anchor", "middle")
            .text( 'Number of Alleles')

        @_getCounts(searchResults)

    _getCounts: (searchResults) ->
        values = []
        for genome in Object.keys(searchResults)
            for gene in Object.keys(searchResults[genome])
                values.push(searchResults[genome][gene])

        console.log("values", values)
        @updateHistogram(values)


    updateHistogram: (values) ->
        histData = @histogram(values)

        # Max y values are increased in discrete steps so that
        # scale is more consistent between updates
        steps = [10,50,100,200,500,800,1000,1200,1500,2000,5000,8000,10000,20000,50000,80000,100000]
        maxSteps = steps.length
        maxY = d3.max(histData, (d) -> d.y)
        yTop = NaN
        for i in [0..maxSteps] by 1
            if maxY < steps[i]
                yTop = steps[i]
                break

        @y = d3.scale.linear()
            .domain([0, yTop])
            .range([@height, 0])

        svgBars = @canvas.selectAll("g.histobar")
            .data(histData)

        # Update existing
        svgBars.attr("transform", (d) => "translate(" + @x(d.x) + "," + @y(d.y) + ")" )

        svgBars.select("rect")
            .attr("x", 0)
            .attr("width", @x.rangeBand())
            .attr("height", (d) => @height - @y(d.y) )
          
        svgBars.select("text")
            .attr("dy", ".75em")
            .attr("y", -14)
            .attr("x", @x.rangeBand() / 2)
            .attr("text-anchor", "middle")
            .text( (d) =>
                if d.y > 0
                    @formatCount(d.y)
                else
                    '' )
            
        # Remove old
        svgBars.exit().remove()
        
        # Insert new
        newBars = svgBars.enter().append("g")
            .attr("class", "histobar")
            .attr("transform", (d) => "translate(" + @x(d.x) + "," + @y(d.y) + ")" )
        
        newBars.append("rect")
            .attr("x", 0)
            .attr("width", @x.rangeBand())
            .attr("height", (d) => @height - @y(d.y) )
            
        newBars.append("text")
            .attr("dy", ".75em")
            .attr("y", -14)
            .attr("x", @x.rangeBand() / 2)
            .attr("text-anchor", "middle")
            .text( (d) =>
                if d.y > 0
                    @formatCount(d.y)
                else
                    '' )
              
        true

