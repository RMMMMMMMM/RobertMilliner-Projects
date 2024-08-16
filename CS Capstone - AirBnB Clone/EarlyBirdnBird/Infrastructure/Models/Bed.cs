using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace Infrastructure.Models
{
    public class Bed
    {
        [Key]
        public int Id { get; set; }

        [Required]
        public string?BedType { get; set; }


        public ICollection<RentalBed> RentalBeds { get; set; }
    }
}
