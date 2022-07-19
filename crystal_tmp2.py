import Crystal as Cr

mc=Cr.Crystal()
mc.load('/dls_sw/i16/software/python/crystal/test.json')
mc.reflection_list(5,refl='sym',anomalous_flag=False,sort='two_theta',print_list=True)

mc.reflection_list(5)

ll=mc.reflection_list(8)

print [list(r[0]) for r in mc.reflection_list(8)]


#doc string
#list of dictionaries?

