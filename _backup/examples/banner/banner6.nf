

params.message = "hello"

data_ch = Channel.from(params.message)

process banner {
  input:
    val(message) from data_ch
  output:
    file(result) into banner_ch
  script:
    result = "banner.txt"
    """
    banner ${params.message} >  $result
    """
}


process countRows {
  validExitStatus 0, 1
  input:
    file(result) from banner_ch
  each col from 1..80
  output:
    stdout col_ch
  script:
  """
    cut -b $col $result > col.txt
    grep -c \\#  col.txt
  """
}

process combineResults {
  input:
  val(res) from col_ch.toList()
  output:
    set  val(size), file("result") into comb_ch
    val(res) into another_ch
  script:
    size = res.size()
    data = res*.trim().join(" ")
    """
     bug
     count_num.py $data > result
   """
}
comb_ch.subscribe { println it }
