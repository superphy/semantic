#Sortable Table
#http://lhorie.github.io/mithril-blog/vanilla-table-sorting.html

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

class Sorting_Table

	controller: () ->
		maxVal = 10
		stepSize = 1
		list = []
		for scale in [maxVal..0] by -stepSize
			console.log scale
			list[scale] = {name:scale, age: maxVal - scale, three: "Hello" + scale}
		{list : list}

	view: (controller) ->
		headers = ['name','age','three']
		m 'table', sort_table(controller.list),
		[
			m 'tr', [
				for x in headers
					m 'th[data-sort-by=' + x + ']', x
			]
			for person in controller.list
				m 'tr', [
					for x in headers
						m 'td', person[x]
				]
		]
